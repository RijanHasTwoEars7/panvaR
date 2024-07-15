#!/bin/bash

# check if the dependencies are present and available on path
command -v tabix >/dev/null || { echo "This shell does not have access to tabix. This software requires tabix in path"; exit; }
command -v vcftools >/dev/null || { echo "This shell does not have access to vcftools. This software requires vcftools in path"; exit; }

# defining defaults
distance=500000 # this will be used if the user does not supply any defaults
window=500000
output_file="panvar_run.txt"
r2_threshold=0.1
minor_allele_frequency=0.05
# defining functions

gene_region_finder () {
  # re-assign the variables here to use inside the function
  chromosome="$1"
  vcf_file="$2"
  loci="$3"
  output_file="$4"
  distance="$5"
  window="$6"
  r2_threshold="$7"
  minor="$8"

  # Check if the output file's directory exists, create if not
  output_dir=$(dirname "$output_file")
  if [[ ! -d "$output_dir" ]]; then
    mkdir -p "$output_dir"
  fi

  # get file base_names for later use
  base_name="$(basename -- $vcf_file)"
  base_name="${base_name%.*}" # this should be the base of the file without the extension

  # Check if the .tbi file exists
  if [[ ! -f "${vcf_file}.tbi" ]]; then
    # Check if a bulk file is supplied by checking the input_file variable
    if [[ -n "$input_file" ]]; then
      echo "Bulk file detected. Generating .tbi file for $vcf_file."
      tabix -p vcf $vcf_file
    else
      read -p "The .tbi file for the given vcf file does not exist. Would you like to generate it? (y/n) " yn
      case $yn in
          [Yy]* ) tabix -p vcf $vcf_file;;
          [Nn]* ) exit;;
          * ) echo "Please answer yes or no."; exit;;
      esac
    fi
  fi


  # Crunching the numbers for the linkage loci range from the input
  snp_start_ld=$((loci-distance))
  snp_stop_ld=$((loci+distance))

  ## make sure that the start LD is not below zero

  if (( ${snp_start_ld} < 1)); then
    snp_start_ld=0
    echo "the given loci was too close to the start and needed to be reset to 0"
  fi

  # generate a base name for the output file using the base name of the input file and the ld range
  tabix_output_file_name="${base_name}_${snp_start_ld}_${snp_stop_ld}"

  tabix_output_file="${output_dir}/${tabix_output_file_name}.txt"

  tabix -h ${vcf_file} ${chromosome}:${snp_start_ld}-${snp_stop_ld} > ${tabix_output_file}

  cat <<- HERE > "${chromosome}_${loci}_temp_region_file.txt"
  Chr Snp
  $chromosome $loci
HERE

  vcftools --vcf $tabix_output_file \
  	--maf $minor \
    --geno-r2-positions "${chromosome}_${loci}_temp_region_file.txt" \
    --ld-window-bp $window \
    --out $output_dir/${tabix_output_file_name}
  
  # remove unwanted files here
  rm "${chromosome}_${loci}_temp_region_file.txt" # the temp file that holds the `--geno-r2-positions` input data for vcftools, was created using the HERE doc.
  
  rm $output_dir/${tabix_output_file_name}.log # the log file that is completely unnecessary for us.

  rm $tabix_output_file

  vcf_file="$output_dir/${tabix_output_file_name}.list.geno.ld"
  
  awk_file_output="$output_dir/${base_name}_${snp_start_ld}_${snp_stop_ld}.ld_data"

  # The sorting should done by the column that has the values for BP, in this case that is the fourth COlumn hence `k4,4n`
  awk -v r2_threshold=$r2_threshold '$6+0 >= r2_threshold'  $vcf_file | grep -v "nan" | grep -v "N_INDV" | sort -k4,4n > $awk_file_output

  # Print the first item of the 4th column
  start=$(awk '{if(NR==1) print $4}' $awk_file_output)
  
  # Print the last item of the 4th column
  stop=$(awk '{last=$4} END{print last}' $awk_file_output)

  ## Format output

  # Check if the output file exists and has the header line
  if [[ ! -f "$output_file" ]]; then
    # If the file doesn't exist, create it and add the header with tabs
    printf "VCF_file\tChrom\tlocus\tdistance\tstart\tstop\tr2_threshold\tminor_allele_freq\n" > "$output_file"
  else
    # If the file exists, check for the header
    if ! grep -q $'VCF_file\tChrom\tlocus\tdistance\tstart\tstop\tr2_threshold\tminor_allele_freq\n' "$output_file"; then
      # If the header is not found, add it with tabs
      sed -i '1iVCF_file\tChrom\tlocus\tdistance\tstart\tstop\tr2_threshold\tminor_allele_freq\n' "$output_file"
    fi
  fi

  # This is where you will write the output of the current iterations
  # The `%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n` formats the outputs to have proper tab spacing
  # Each item should have its wown tab wiht and \t and I believe the %s is for type
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$base_name" "$chromosome" "$loci" "$distance" "$start" "$stop" "$r2_threshold" "$minor" >> "$output_file"
  
  rm $vcf_file # remove unnecessary intermediate file
  rm $awk_file_output
}

# CLI flags
while test $# -gt 0; do
  case "$1" in
    -f|--file)
      shift
      input_file="$1"
      shift
      ;;
    -o|--output)
      shift
      output_file=$1
      shift
      ;;
    *)
      # Handle other flags or break out of the loop
      break
      ;;
  esac
done

# Process the input file or command line arguments
if [[ -n "$input_file" ]]; then
  if [[ ! -f "$input_file" ]]; then
    echo "Input file '$input_file' not found."
    exit 1
  fi

  # Process the SSV file
while IFS=$' ' read -r chromosome vcf_file loci custom_distance custom_window custom_r2_threshold custom_minor_allele_frequency; do # this is where the file reads things
  # Skip the loop iteration if the line is empty
  [[ -z "$chromosome" ]] && continue
  # Use the provided distance and window if present, otherwise use the defaults
  gene_region_finder "$chromosome" "$vcf_file" "$loci" "$output_file" "${custom_distance:-$distance}" "${custom_window:-$window}" "${custom_r2_threshold:-$r2_threshold}" "${custom_minor_allele_frequency:-$minor_allele_frequency}"
done < "$input_file"


else
  # Reset the argument pointer if necessary
  # set -- 

  while test $# -gt 0; do
    case "$1" in
      -c|--chromosome)
        shift
        chromosome=$1
        shift
        ;;
      -v|--vcf_file)
        shift
        vcf_file=$1
        shift
        ;;
      -o|--output)
        shift
        output_file=$1
        shift
        ;;
      -l|--loci)
        shift
        loci=$1
        shift
        ;;
      -d|--distance)
        shift
        distance=$1
        shift
        ;;
      -w|--window)
        shift
        window=$1
        shift
        ;;
      -r2|--r2_threshold)
        shift
        r2_threshold=$1
        shift
        ;;
      -m|--minor)
        shift
        minor_allele_frequency=$1
        shift
        ;;
      *)
        # Handle unknown option or break
        break
        ;;
    esac
  done
  gene_region_finder "$chromosome" "$vcf_file" "$loci" "$output_file" "$distance" "$window" "$r2_threshold" "$minor_allele_frequency"
fi
