#!/bin/bash

# Check for dependencies
if ! command -v vcfEffOnePerLine.pl &> /dev/null
then
    echo "vcfEffOnePerLine.pl could not be found"
    exit 1 # Changed exit code to 1 for error indication
fi

if [ ! -f "./snpSift.jar" ]; then
    echo "snpSift.jar does not exist in the current directory"
    exit 1 # Changed exit code to 1 for error indication
fi

process_gene(){

    gene_name="$1"
    vcf_file="$2"
    output_file="$3"
    gene_location="$4"
    
    # Check if the output file's directory exists, create if not
    output_dir=$(dirname "$output_file")
    if [[ ! -d "$output_dir" ]]; then
      mkdir -p "$output_dir"
    fi

    # Removed the commented-out section for creating a sub-directory for the gene

    # Check if gene_name exists in gene_location file
    gene_line=$(grep -F "$gene_name" "$gene_location")
    if [[ -z "$gene_line" ]]; then
        echo "Warning: $gene_name does not exist in $gene_location"
        return 1
    fi

    # If gene_name exists in gene_location file, use awk to create the required string
    gene_string=$(echo "$gene_line" | awk '{print $1 ":" $2 "-" $3}') # this creates the query for tabix to use against the vcf file

    # Corrected the use of variables in tabix command
    tabix -h "$vcf_file" "$gene_string" | \
        grep -E "#|${gene_name}" > "${gene_name}_region_filtered.vcf"

    # Corrected the use of variables in tabix command
    tabix -H "$vcf_file" | tail -n 1 > "${gene_name}.vcf_header"

    # New line added here
    cat "${gene_name}_region_filtered.vcf" | vcfEffOnePerLine.pl | grep -E "#|${gene_name}" | \
    java -jar snpSift.jar extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" | \
    grep 'Chr' >> "${gene_name}_processed.vcf"

    # Rijan: remove multiple/non-relevant gene snps
    grep -v "-" "${gene_name}_processed.vcf" > "${gene_name}_single_processed.vcf"

    # Improved sed command readability
    sed -e 's|0/0|0|g' \
        -e 's|0/1|1|g' \
        -e 's|1/1|2|g' \
        -e 's|./.|NA|g' "${gene_name}_single_processed.vcf" > "${gene_name}_numerical_allele_scores.vcf"

    # set up the output file if need be
    header_file="${gene_name}.vcf_header"
    impactful_file="$output_file.impactful_hits.vcf"

    # Check if the impactful hits file exists before reading
    if [[ -f "$impactful_file" ]]; then
        read -r first_line < "$impactful_file"
    else
        touch "$impactful_file" # Create the file if it doesn't exist
        first_line=""
    fi

    # Read the content of the header file
    read -r header < "$header_file"

    # Check if the first line matches the header
    if [ "$first_line" != "$header" ]; then
      # Add the header to the top of the impactful hits file
      { echo "$header"; cat "$impactful_file"; } > temp_file && mv temp_file "$impactful_file"
    fi

    # append outcome to file
    grep -E 'CHROM|LOW|MODERATE|HIGH' "${gene_name}_numerical_allele_scores.vcf" >> "$output_file"

    # remove intermediate files
    rm "${gene_name}.vcf_header"
    rm "${gene_name}_region_filtered.vcf"
    rm "${gene_name}_single_processed.vcf"
    rm "${gene_name}_numerical_allele_scores.vcf"
    rm "${gene_name}_processed.vcf"
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
      output_file="${1:-"impactful_hits"}"
      shift
      ;;
    -loc|--gene_location) # This is the file with the gene name locations extended
      shift
      GENE_LOCATION="$1"
      shift
      ;;
    -v|--vcf) # This is the vcf file that needs to be searched
      shift
      VCF_FILE="$1"
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
  while IFS=$' ' read -r gene_name custom_gene_location_file custom_vcf_file custom_output; do
    # Skip the loop iteration if the line is empty
    [[ -z "$gene_name" ]] && continue # Corrected the variable name from chromosome to gene_name
    # Use the provided distance and window if present, otherwise use the defaults
    # position and order matters when calling a function so be careful
    process_gene "$gene_name" "${custom_vcf_file:-$VCF_FILE}" "${custom_output:-$output_file}" "${custom_gene_location_file:-$GENE_LOCATION}"  
  done < "$input_file"
else

  # Reset the argument pointer if necessary
  # set -- 

  while test $# -gt 0; do
      case "$1" in
        -g|--gene_name) # This is the name of the gene that you want searched
          shift
          GENE_NAME="$1"
          shift
          ;;
        -loc|--gene_location) # This is the file with the gene name locations extended
          shift
          GENE_LOCATION="$1"
          shift
          ;;
        -v|--vcf) # This is the vcf file that needs to be searched
          shift
          VCF_FILE="$1"
          shift
          ;;
        -o|--output) # This is the directory where the output will be piped to
          shift
          output_file="$1"
          shift
          ;;
        *)
          # Handle unknown option or break
          break
          ;;
      esac
    done
  # Check if the gene_location file exists
  if [[ ! -f "$GENE_LOCATION" ]]; then
      echo "The gene_location file does not exist. Please enter the correct path:"
      read GENE_LOCATION
  fi

  # Check if the vcf_file exists
  if [[ ! -f "$VCF_FILE" ]]; then
      echo "The vcf file does not exist. Please enter the correct path:"
      read VCF_FILE
  fi

  # Check if the output_dir exists
  output_dir=$(dirname "$output_file")
  if [[ ! -d "$output_dir" ]]; then
      echo "The output directory does not exist. Do you want to create it? (yes/no):"
      read RESPONSE
      if [[ "$RESPONSE" == "yes" ]]; then
          mkdir -p "$output_dir"
          echo "Output directory created at $output_dir"
      else
          echo "Exiting the program."
          exit 1
      fi
  fi

  process_gene "$GENE_NAME" "$VCF_FILE" "$output_file" "$GENE_LOCATION"  

fi