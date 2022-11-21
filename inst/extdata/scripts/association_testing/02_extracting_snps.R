#! /usr/bin/env Rscript
'
This script uses bgenix or PLINK 2 to extract SNPs from genetic data, in .bgen, .bed or .pgen formats, using SNP identifiers and chromosome number. As such, bgenix and/or PLINK should be installed.
First, the genetic_file_guide_template.csv file must be edited to provide the file location of the genetic files and corresponding .sample/.fam/.psam file for each of the 22 autosomal and 2 sex chromosomes. If genetic data for a particular chromosome or chromosome(s) are missing, the corresponding rows in the data/genetic_file_guide_template.csv file should be deleted. If the genetic data are not stored per chromosome, then different values can be used in the chromosome column. The genetic_file_guide_template.csv is found within the R package in the extdata folder compressed as a gzip file. A template for SNP_list.csv is also saved in extdata.

Usage:
    02_extracting_snps.R (--genetic_file_guide=<FILE> --SNP_list=<FILE> --analysis_folder=<FOLDER>) (--plink_input | --bgen_input) [--plink_exe=<text> --bgenix_exe=<text> --plink_type=<text> --ref_bgen=<text> --variant_save_name=<name> --no_delete_temp]

Options:
    -h --help  Show this screen.
    -v --version  Show version.

    Mandatory inputs
    --genetic_file_guide=<FILE>           Full path of the completed genetic_file_guide_template.csv.

    --SNP_list=<FILE>                     Full path of the file describing which SNPs are to be extracted from the genetic data files ahead of association testing.

    --analysis_folder=<FOLDER>            Full path of the folder that will contain the data for the SNPs given by SNP_list. A temporary folder, named "temp_plink"
                                          by default, will be created within the folder specified by this argument as a place to hold temporary files needed for the
                                          SNP extraction process. Unless otherwise requested by specifying the no_delete_temp argument, the temporary folder and its
                                          contents will be deleted following successful SNP extraction.
    Select one of
    --bgen_input                          Specify that the genetic data files are in .bgen format.
    --plink_input                         Specify that the genetic data files are in PLINK format (.bed or .pgen).

    Options
    --plink_exe=<text>                    Full path to the PLINK2 executable. [default: plink2]
    --plink_type=<text>                   Specify whether the PLINK-formatted genetic data are in .bed or .pgen format. [default: bed]
    --ref_bgen=<text>                     One of the following three values that specifies which allele is to be used as the reference: ref-first (first allele is the reference, default),
                                          ref-last (last allele is the reference), red-unknown (last allele is provisionally treated as the reference). [default: ref-first]
    --bgenix_exe=<text>                   Full path to the bgenix executable. [default: bgenix]
    --variant_save_name=<name>            Name of the output genetic data files. [default: variants_for_association]
    --no_delete_temp                      Specify whether the temporary folder containing intermediate files should be retained (TRUE) or deleted (FALSE). Default is TRUE.

' -> doc


suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.2 02_extracting_SNPs')

library(DeepPheWAS)

snp_extractor(genetic_file_guide=arguments$genetic_file_guide,
              SNP_list=arguments$SNP_list,
              analysis_folder=arguments$analysis_folder,
              bgen_input=arguments$bgen_input,
              plink_input=arguments$plink_input,
              plink_exe=arguments$plink_exe,
              plink_type=arguments$plink_type,
              ref_bgen=arguments$ref_bgen,
              bgenix_exe=arguments$bgenix_exe,
              variant_save_name=arguments$variant_save_name,
              no_delete_temp=arguments$no_delete_temp)
