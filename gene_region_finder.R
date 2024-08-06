#!/bin/Rscript
options(scipen=999) # This makes sure that shell commands do not get converted to scientific notation.

# Dependencies
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

# Define the command-line arguments
option_list <- list(
  make_option(c("-d", "--distance"), type="double", default=500000, help="The upstream/distance from the loci for which the VCF will be subsetted?"),
  
  make_option(c("-w", "--window"), type="double", default=500000, help="What is the LD window for LD calculation?"),
  
  make_option(c("-r", "--r2_threshold"), type="double", default=0.1, help="What is the r^2 threshold for the locus pair?"),
  
  make_option(c("-o", "--output_file"), type="character", default="panvar_run.txt", help="Where should the output file be stored?"),
  
  make_option(c("-c", "--chromosome"), type="character", help="The chromosome in the file of interest. Typos matter!"),

  make_option(c("-v", "--vcf_file"), type="character", help="The path to the VCF file of interest."),

  make_option(c("-l", "--loci"), type="integer", help="The loci of interest. The first locus is used as the reference."),
	
	make_option(c("-b","--bulf_file"), type = "character", default = NULL,help = "A CSV or TSV file that has bulk inputs.")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Functions
# Check if external bin is on path
is_bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

# define a function that take a path to a vcf file and checks if the .tbi file exists,
# if it exists return TRUE, else
# If not ask the user if they want to generate a .tbi file and if the user says no quit the script

proper_tbi <- function(vcf_file) {
  if (!file.exists(paste0(vcf_file, ".tbi"))) {
    # ask the user if they want to generate the .tbi file
    cat("The .tbi file for the given vcf file does not exist. Would you like to generate it? (y/n) ")
    answer <- tolower(readLines("stdin", n = 1))
    if (answer == "y") {
      # generate the .tbi file
      print("Asking tabix to generate the .tbi file...")
      system(paste0("tabix -p vcf ", vcf_file))
      return(TRUE)
    } else {
      stop("The .tbi file for the given vcf file does not have a tabix index. Please generate it and try again.")
    }
  } else {
    return(TRUE)
  }
}
# We need to check if bedtools and tabix are installed
# To check if a binary exists in R use this, from: https://search.r-project.org/CRAN/refmans/bedr/html/check.binary.html

if (!is_bin_on_path("tabix")) {
  stop("Tabix is not accessible to this shell. Please install it or make it accessible and try again.")
}

if (!is_bin_on_path("vcftools")) {
  stop("bedtools is not accessible to this shell. Please install it and try again.")
}

gene_region_finder <- function(chromosome, vcf_file, loci, output_file, distance, window, r2_threshold) {


	if (!file.exists(vcf_file)){
		stop("The VCF file is not accessible from this location. Did you supply the whole path?")
	}
  # Check if the output file's directory exists, create if not
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # get file base_names for later use
  base_name <- basename(vcf_file)
  base_name <- sub("\\..*$", "", base_name) # this should be the base of the file without the extension

  # check if the .tbi file exists for the current vcf file and generate it if not
  proper_tbi(vcf_file)

  # Crunching the numbers for the linkage loci range from input
  snp_start_ld <- loci - distance
  snp_end_ld <- loci + distance

  ## make sure that the start loci is not less than 1 
  if (snp_start_ld < 1) {
    snp_start_ld <- 0
    print("The given loci generated a start point less than 1. Setting it to 0.")
  }

  # generate the base name for the output file using the base name of the input file and ld range
  tabix_file_subset <- paste0(base_name, "_", snp_start_ld, "_", snp_end_ld,".txt")

  # asking tabix to subset the file
  system(paste0("tabix -h ", vcf_file, " ", chromosome, ":", snp_start_ld, "-", snp_end_ld, " > ", tabix_file_subset))

  # check to see if the tabix output file has at least 1 line
  # that way we can know that the subsetting did not go awray
  if (length(readLines(tabix_file_subset)) == 0) {
    stop("The program was unable to subset your VCF file for the given input. Please check that the right input file, loci and chromsome name were supplied. Checking for typos might be a good idea, Chr_05 is not the same as Chr05.")
  }

  # Here doc file for tabix snp calculation
  dt <- data.table(
    Chr = c(chromosome),
    Snp = c(loci)
  )

  dt %>%
    fwrite("current_geno_r2_positions_table.txt", sep = "\t", col.names = TRUE)

  # Make a call to the right file

  system2("echo",args =c(window))

  system2("vcftools",
  args = c(
    "--vcf", tabix_file_subset,
    "--geno-r2-positions", "current_geno_r2_positions_table.txt", # Because the file name is constant we do not need a variable here
    "--ld-window-bp", window,
    "--out", file.path(output_dir, tabix_file_subset)
    )
  )

  current_ld_file_path = paste0(output_dir, "/", tabix_file_subset, ".list.geno.ld")

  # check the output of vcftools for length to see that it isn't empty
  if (length(readLines(current_ld_file_path)) == 0) {
    stop("The program was unable to calculate LD for the given loci. Please check why this happened. The program needs LD to generate the final reports.")
  }

  # TODO: delete log file and maybe move subset file
  # read the output of LD calulation into a data table

  current_ld_calculations <- fread(current_ld_file_path, header = TRUE)

  # filter and sort
  sorted_by_ld <- current_ld_calculations %>%
  filter(`R^2` >= r2_threshold, 
         !is.na(`R^2`), 
         `R^2` != "N_INDV") %>%
  arrange(`R^2`)
  
  start <- first(sorted_by_ld$POS2)

  stop <- last(sorted_by_ld$POS2)

  # add current data in a table format

  # assuming new_data is a data.table with the same structure as current_table
  current_input <- data.table(
    VCF_file = c(vcf_file),
    Chrom = c(chromosome),
    locus = c(loci),
    distance = c(distance),
    r2_threshold = c(r2_threshold),
    start = c(start),
    stop = c(stop)
  )

  # check if the output file exists and
  # if it doesn't create a new table
  # if it does read it in
  if (file.exists(output_file)) {
    current_table <- fread(output_file, header = TRUE)

    current_table <- rbind(current_table, current_input)

    deduplicated_table <- current_table[!duplicated(current_table),]

    deduplicated_table %>% 
      fwrite(output_file, sep = "\t", col.names = TRUE)
  } else{
  
    current_table <- current_input

    current_table %>% 
      fwrite(output_file, sep = "\t", col.names = TRUE)
  }
}

# calling the function for a sample run

if (opt$bulk_file == NULL) {
	
	gene_region_finder(
	chromosome = opt$chromosome,
	vcf_file = opt$vcf_file,
	loci = opt$loci,
	output_file = opt$output_file,
	distance = opt$distance,
	window = opt$window,
	r2_threshold = opt$r2_threshold 
	)
	
} else {

	bulk_table = fread(opt$bulk_file)

	right_columns = c("VCF_file","Chromosome","Loci","Distance","Window","R2_threshold","Output_file") # These are the right columns
	
	if (nrows(bulk_table) <= 0){
		stop("This table is either an empty file or has no values past metadata. Please use --sample_table to see what a proper sample table looks like. Spelling and capitals matter.")
	}

	current_columns = colnames(bulk_table)

	if(!all(current_columns %in% right_columns)){
		
		print("The bulk table does not have all the right columns. They are either absent or mis-spelled. The following are the columns you should have.")

		print(right_columns)

		stop()
	}
	
	for (i in 1:nrow(bulk_table){
		# I don't like that I have to spell all columns out by hand but this is the best I know for now.
		# Now mixing and matching things is odd because that would need try catches
		gene_region_finder(
			chromosome = bulk_table[i,]$Chromosome,
			vcf_file = bulk_table[i,]$VCF_file,
			loci = bulk_table[i,]$Loci,
			distance = bulk_table[i,]$Distance,
			window = bulk_table[i,]$Window,
			r2_threshold = bulk_table[i,]$R2_threshold,
			output_file = opt$Output_file
		)
	}
}		 