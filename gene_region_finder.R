#!/bin/Rscript

# Dependencies
library(data.table)
library(tidyverse)
library(optparse)

# Define the command-line arguments
option_list <- list(
  make_option(c("-d", "--distance"), type="double", default=500000, help="The upstream/distance from the loci for which the VCF will be subsetted?"),
  
  make_option(c("-w", "--window"), type="double", default=500000, help="What is the LD window for LD calculation?"),
  
  make_option(c("-r", "--r2_threshold"), type="double", default=0.1, help="What is the r^2 threshold for the locus pair?"),
  
  make_option(c("-o", "--output_file"), type="character", default="panvar_run.txt", help="Where should the output file be stored?")
  
  make_option("-c", "--chromosome", type="character", help="The chromosome in the file of interest. Typos matter!")

  make_option("-v", "--vcf_file", type="character", help="The path to the VCF file of interest.")

  make_option("-l", "--loci", type="integer", help="The loci of interest. The first locus is used as the reference.")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Calculate the sum
sum <- opt$num1 + opt$num2

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
    answer <- readline(paste0("The .tbi file for the given vcf file does not exist. Would you like to generate it? (y/n) "))
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
  stop("Tabox is not accessible to this shell. Please install it and try again.")
}

if (!is_bin_on_path("vcftools")) {
  stop("bedtools is not accessible to this shell. Please install it and try again.")
}

gene_region_finder <- function(chromosome, vcf_file, loci, output_file, distance, window, r2_threshold) {

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
  tabix_output_file_name <- paste0(base_name, "_", snp_start_ld, "_", snp_end_ld,".txt")

  # asking tabix to subset the file
  system(paste0("tabix -h ", vcf_file, " ", chromosome, ":", snp_start_ld, "-", snp_end_ld, " > ", tabix_output_file_name))

  # check to see if the tabix output file has at least 1 line
  # that way we can know that the subsetting did not go awray
  if (length(readLines(tabix_output_file_name)) == 0) {
    stop("The program was unable to subset your VCF file for the given input. Please check that the right input file, loci and chromsome name were supplied. Checking for typos might be a good idea, Chr_05 is not the same as Chr05.")
  }

  # Here doc file for tabix snp calculation
  chromosome <- 
  dt <- data.table(
    Chr = c(chromosome),
    Snp = c(loci)
  )

  dt %>%
    fwrite("current_geno_r2_positions_table.txt", sep = "\t", col.names = TRUE)

  # Make a call to the right file

  system2("vcftools",
  args = c(
    "--vcf", tabix_output_file,
    "--geno-r2-positions", "current_geno_r2_positions_table.txt", # Because the file name is constant we do not need a variable here
    "--ld-window-bp", window,
    "--out", file.path(output_dir, tabix_output_file_name)
    )
  )

  current_ld_file_path = paste0(output_dir, "/", tabix_output_file_name, ".list.geno.ld")

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
  
  start <- sorted_by_ld$POS2[1]

  stop <- last(sorted_by_ld$POS2)

  # add current data in a table format

  # assuming new_data is a data.table with the same structure as current_table
  current_input <- data.table(
    VCF_file = c(new_VCF_file),
    Chrom = c(new_chromosome),
    locus = c(new_loci),
    distance = c(new_distance),
    r2_threshold = c(new_r2_threshold),
    start = c(new_start),
    stop = c(new_stop)
  )

  # check if the output file exists and
  # if it doesn't create a new table
  # if it does read it in
  if (file.exists(output_file)) {
    current_table <- fread(output_file, header = TRUE)

    current_table <- rbind(current_table, current_input)

    current_table %>%
      fwrite(output_file, sep = "\t", col.names = TRUE)
  } else{
    current_table <- current_input

    current_table %>%
      fwrite(output_file, sep = "\t", col.names = TRUE)
  }
}