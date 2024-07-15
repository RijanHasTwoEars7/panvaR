
Disclaimer: This README is a WIP.

# Getting started

You can clone this repository with the commands `git clone https://gitlab.com/ddpsc/popgen/panvar.git`. <- This command will clone the main branch of `Panvar`, which is the branch meant for the end user.

## Using panvar

Panvar is divided into four steps:

1. Window determination
2. List genes
3. List polymorphisms
4. Phenotype statistics

### Window determination

The task of window determincation is done by the script, creatively named, `gene_region_finder.sh`. This bash script is designed to identify genomic regions around a specified locus within a VCF (Variant Call Format) file using the `tabix` and `vcftools` utilities. It allows users to specify a chromosome, a locus, and a distance to define the region of interest.

**Dependencies**
- `tabix`: Used for indexing and retrieving parts of the VCF file.
- `vcftools`: Used for calculating linkage disequilibrium (LD) and filtering data.

**Installation**
Ensure that both `tabix` and `vcftools` are installed and available in your system's PATH.

**Usage**
The script can be run with either command-line arguments or an input file containing the parameters.

**Command-Line Arguments**
- `-c` or `--chromosome`: Specify the chromosome number.
- `-v` or `--vcf_file`: Specify the path to the VCF file.
- `-l` or `--loci`: Specify the locus position.
- `-d` or `--distance`: Specify the distance from the locus to define the region (default: 500000).
- `-w` or `--window`: Specify the window size for LD calculation (default: 500000).
- `-o` or `--output`: Specify the output file name (default: "panvar_run.txt").
- `-r2`: Specify the r2 cut-off for LD scores.

**Input File Mode**
- `-f` or `--file`: Specify the path to the input file containing the parameters.
- The input file should be a space-separated values (SSV) file with the following columns, in this order: chromosome, VCF file path, locus, distance, and window.
- The input file's rows can skip distance and window, in which case the defaults would be used.

Here is what an actual input file can look like:

(**Note**:Header only included for Markdown formatting purposes.)

| Chr | path_to_vcf_file | 59629117 | 500000 | 500000 |
|:-----:|:----------------:|:--------:|:------:|:------:|
| Chr06 | path_to_vcf_file | 59629117 | 500000 | 500000 |
| Chr08 | path_to_vcf_file | 308199   | 500000 | 500000 |
| Chr06 | path_to_vcf_file | 59629117 | 500000 | 500000 |
| Chr06 | path_to_vcf_file | 59629117 | 500000 | 500000 |
| Chr08 | path_to_vcf_file | 308199   | 500000 | 500000 |
| Chr02 | path_to_vcf_file | 4811079  | 500000 | 500000 |

**Output**
The script generates several files during its operation, including:
- A `.tbi` index file for the VCF file (if not present).
- A temporary region file used for LD calculations.
- An output file containing the LD data and the defined genomic region.

Here is what a typical output from `gene_region_finder.sh` will look like:

| VCF_file                                 | Chrom | locus    | distance | start    | stop     |
|------------------------------------------|-------|----------|----------|----------|----------|
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr08.snp_final.vcf  | Chr08 | 308199   | 500000   | 17553    | 804806   |
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr08.snp_final.vcf  | Chr08 | 308199   | 500000   | 17553    | 804806   |
| Sorghum_V5.maf0.001.Chr02.snp_final.vcf  | Chr02 | 4811079  | 500000   | 4312510  | 5299815  |

**Examples**
To run the script with command-line arguments:
```bash
./gene_region_finder.sh -c 1 -v sample.vcf -l 123456 -d 100000 -w 100000 -o output.txt
```

To run the script with an input file:
```bash
./gene_region_finder.sh -f input_parameters.ssv -o output.txt
```

**Note**: The script checks for the presence of the `.tbi` index file and prompts the user to create one if it's missing. It also creates necessary directories for the output files.

**Troubleshooting**
- Ensure `tabix` and `vcftools` are correctly installed and accessible.
- Verify the format and existence of the input VCF file and the `.tbi` index file. Depending on how your `tabix` and `vcftools` were compiled, you might suffer buffer overflows. 
- Check the permissions of the directory where the output files will be saved.

### Listing genes

**User Manual for Gene Region Parser Script**

**Introduction**
This Python script is designed to parse gene location data and extract specific gene regions based on user-defined chromosome coordinates. It can process individual regions or bulk input from a TSV file.

**Requirements**
- Python 3.x
- Pandas library

**Installation**
Ensure Python and Pandas are installed on your system. If Pandas is not installed, you can install it using pip:
```bash
pip install pandas
```

**Usage**
The script can be run from the command line with the following arguments:

- `-l` or `--gene_location_file`: Path to the gene location file (required).
- `-start` or `--start_point`: Start point of the gene region (integer, required).
- `-stop` or `--stop_point`: Stop point of the gene region (integer, required).
- `-c` or `--chromosome`: Chromosome of the gene region (required).
- `-o` or `--output_file`: Path to the output file where results will be saved (required).
- `-b` or `--bulk_input_file`: Path to the bulk input file as a TSV (optional).

**Example Command**
```bash
python script.py -l gene_locations.tsv -start 100000 -stop 200000 -c chr1 -o output_directory/output_file
```

**Bulk Processing**
To process multiple regions at once, provide a TSV file with columns `Chrom`, `start`, `stop`, and optionally `genloc` and `output_file` for custom file paths.

here is what a sample bulk input file looks like:

| VCF_file                                 | Chrom | locus    | distance | start    | stop     |
|------------------------------------------|-------|----------|----------|----------|----------|
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr08.snp_final.vcf  | Chr08 | 308199   | 500000   | 17553    | 804806   |
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr06.snp_final.vcf  | Chr06 | 59629117 | 500000   | 59129242 | 60126036 |
| Sorghum_V5.maf0.001.Chr08.snp_final.vcf  | Chr08 | 308199   | 500000   | 17553    | 804806   |
| Sorghum_V5.maf0.001.Chr02.snp_final.vcf  | Chr02 | 4811079  | 500000   | 4312510  | 5299815  |

Notice that having extraneous fields is fine so long as your supplied table has `Chrom`, `start` and `stop` fields.

You would call the above table as follows:
`python gene_list_producer.py --bulk_input_file <bulk_input_table> --gene_location_file <gene_location_data> --output_file <output_file_path>`

**Output**
The script generates two files:
- A CSV file containing the full gene table.
- A text file listing unique gene names.

**Function Descriptions**
- `region_file_parser`: Reads the gene location file and filters genes within the specified region.
- `write_output_files`: Writes the output files to the specified directory.

**Error Handling**
Ensure all required arguments are provided and that file paths are correct. The script will raise errors if files cannot be found or arguments are missing.

### List polymorphisms

**User Manual for Gene Processing Bash Script**

### Introduction
This script processes genetic data from VCF (Variant Call Format) files. It filters and extracts information based on specific gene names and their locations, then outputs the processed data into a specified file.

### Dependencies
- `vcfEffOnePerLine.pl`: A Perl script to format VCF files into one effect per line.
- `snpSift.jar`: A Java-based tool for filtering and manipulating VCF files.
- `tabix`: A tool that creates an index for a TAB-delimited genome position file and allows fast retrieval of data.

### Installation
Ensure that `vcfEffOnePerLine.pl` is installed and `snpSift.jar` is present in the current directory before running the script.

### Usage
The script can be run for either bulk processing of genes using a tab-separated values (TSV) file or for single gene processing.

#### Bulk Processing
1. Prepare a TSV file with the following columns: gene name, custom VCF file path (optional), custom gene location file path (optional), and custom output file path (optional).
2. Run the script with the `-f` flag followed by the path to the TSV file.

Here is what a typical input file looks like

```
| Sobic.006G235650 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | custom_gene_location (optional) | custom_output (optional) |
| Sobic.006G235700 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G235800 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G235801 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G235900 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G236000 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G236300 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G236400 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G236500 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
| Sobic.006G236600 | /scratch/vcf_files/annotated/sbicolor/Sorghum_V5.maf0.001.Chr06.snp_final.vcf.gz | defaults to `-loc`              | defaults to `-o`         |
```
Notice that the `custom_gene_location` and `custom_output` are optional and default to their flags **when** not in the bulk table's rows.

```bash
bash list_polymorphisms.sh -f <input_table> -loc <gene_location_file> -o <output_name>
```

#### Single Gene Processing
1. Provide the gene name, VCF file path, gene location file path, and output file path using the respective flags: `-g`, `-v`, `-loc`, and `-o`.

```bash
bash list_polymorphisms.sh -g GENE_NAME -v VCF_FILE -loc GENE_LOCATION_FILE -o OUTPUT_FILE
```

### Flags
- `-f` or `--file`: Path to the input TSV file for bulk processing.
- `-g` or `--gene_name`: Name of the gene to process.
- `-v` or `--vcf`: Path to the VCF file.
- `-loc` or `--gene_location`: Path to the gene location file.
- `-o` or `--output`: Path to the output file (default is "impactful_hits").

### Output
The script will output a VCF file with filtered genetic information based on the provided gene names and locations. Intermediate files created during processing are automatically removed.

### Error Handling
The script checks for the existence of required files and dependencies at runtime. If any checks fail, it will output an error message and terminate.

### Notes
- Ensure that the gene location file is formatted correctly with gene names and their corresponding locations.
- The script assumes that the VCF file has an associated `.tbi` index file. If not present, it will attempt to create one using `tabix`.

### Statistical calculations

**User Manual for R Script**

**Introduction**
This R script is designed to perform statistical analysis on genotype and phenotype data. It includes functions for data cleaning, mutation rate calculation, ANOVA testing, and significance adjustment.

**Dependencies**
- `dplyr`: For data manipulation.
- `data.table`: For efficient data handling.
- `argparse`: For parsing command-line arguments.

**Input Requirements**
1. **Genotype File (--geno)**: A file path to the genotype data. The data should be in a tabular format where each row represents an individual and columns represent genotypic information.
2. **Phenotype File (--pheno)**: A file path to the phenotype data. This file should also be in a tabular format with rows representing individuals and at least one column of phenotypic data.
3. **Alpha (--alpha)**: A significance level for statistical tests. The default is set to 0.05.
4. **Output Directory (--output)**: A file path to the desired output directory. If not specified, the current working directory is used.

**Output Description**
The script will generate a summary table with the following columns:
- Metadata fields from the genotype data.
- Adjusted p-values from the ANOVA tests.
- Mutation counts for different genotypes (mut, wild, hybrid).
- A Bonferroni corrected significance value.
- A significance column indicating whether the mutation rates are significant based on the Bonferroni correction.

**Execution Instructions**
To run the script, use the command line to navigate to the directory containing the script and execute it with Rscript, providing the necessary arguments. For example:
```sh
Rscript script_name.R --geno path/to/geno_file --pheno path/to/pheno_file --alpha 0.05 --output path/to/output_directory
```

**Note**: Ensure that the input files exist and are in the correct format before running the script to avoid errors.
