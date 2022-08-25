#!/usr/bin/env Rscript
'

This script is an alternative to the PheWAS_association_PLINK.R script for running association analyses using R. While not streamlined for large-scale PheWAS, this script offers greater flexibility in the kind of analyses that can be run, such as fitting non-standard models or if the genetic data are not in a format readable by standard software. Further, this script allows for analyses of genomic risk scores (GRS), polygenic risk score (PRS), and multiallelic copy number variants (CNV).

GRS/PRS is inputted via a csv file with trait,group,phenotype_data,genetic_data,variable this acts as a map for the functions. See user guide for full description of this input but in brief, it allows the user to input any number of traits that may be analysed in any of the availble groupings (such as ancestry) and point to the location of the appropriate genetic and phenotype files. GRS/PRS data is analysed per-trait-group combination, results are saved as a R list object per trait that contain all trait-group combinations, this allows multiple traits to be analysed across multiple groups in a single input.

Non-GRS/PRS data is assumed to be an alternative genetic measurement with a column of IDs followed by 1-n columns of genetic data, example CNVs. Results are per column of the input file (minus ID column) with a single R object saved containing a list of results per-group such that each non-ID column will be tested for association with all available phenotypes in each of the inputted groups. Unlike the GRS/PRS input the full location of all the phenotype files is required as input.

If a covariate file is included then it must contain a column age, when selecting covariates note currently only participants with full covariate data will be analysed. Other options perform filtering on the phenotypes that are analysed, by default all phenotypes that are held in the respective phenotype tables created by 01_phenotype_preparation.R are analysed. Minimum case numbers are options to allow the input of genetic variables that do not fully overlap with the samples used for making the phenotypes.


Usage:
    03b_R_association_testing.R --analysis_folder=<folder> --GRS_input=<FILE> [--covariates=<FILE> --N_cores=<number> --PheWAS_manifest_overide=<FILE>] [--phenotype_inclusion_file=<FILE> | --phenotype_exclusion_file=<FILE>]

    03b_R_association_testing.R --analysis_folder=<folder> --non_GRS_data=<FILE> --phenotype_files=<FILE> --analysis_name=<name> [--covariates=<FILE> --group_name_overide=<text> --N_cores=<number> --PheWAS_manifest_overide=<FILE>] [--phenotype_inclusion_file=<FILE> | --phenotype_exclusion_file=<FILE>]


Options:
    -h --help  Show this screen.
    -v --version  Show version.

    Mandatory Inputs
    --analysis_folder=<FOLDER>            Full path of the directory in which the association results will be saved.
    One of
    --GRS_input=<FILE>                    Full file path of the GRS_input csv file. See user guide for more information on format.
    --non_GRS_data=<FILE>                 Full file path to the non_GRS genetic data. See user guide for more information on format.

    If non_GRS_data must include
    --phenotype_files=<FILE>              Comma-separated list containing the full paths of the phenotype data.
                                          Only required for non_GRS_data.
    --analysis_name=<name>                Name for the analysis, is used later in saving tables, so should distignuish between other analyses. For
                                          GRS analysis this name is always the trait being analysed.
    Options
    --covariates=<FILE>                   Full file path for the covariates file see user guide for more information.
    --group_name_overide=<text>           Comma-separated list containing alternative group names. By default, group names are extracted from the suffix of
										  the file names provided in the phenotype_files argument. For example, if a file named /home/phenotypes/EUR_phenotypes.csv
										  were provided, the corresponding group name would be "EUR". This argument allows for a different group name to be specified.
										  The order of names provided should match the files specified in the phenotype_files argument.

    --N_cores=<number>                    Number of cores requested if parallel computing is desired. Defaults to single core computing.

    --phenotype_inclusion_file=<FILE>     Full file path to a txt file containing single column containing full PheWAS_ID of phenotypes that will be
                                          included. Cannot be used with phenotype_exclusion_file argument.
    --phenotype_exclusion_file=<FILE>     Full file path to a txt file containing single column containing full PheWAS_ID of phenotypes that will be
                                          excluded. Cannot be used with phenotype_inclusion_file argument.

    --binary_Case_N=<number>              Number that represents the minimum number of cases for binary phenotype inclusion.
                                          [default: 50]
    --quantitative_Case_N=<number>        Number that represents the minimum number of cases for quantitative phenotype inclusion.
                                          [default: 100]
    --PheWAS_manifest_overide=<FILE>      Full file path of the alternative PheWAS_manifest file.

' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.2 03b_R_association_testing.R')

library(DeepPheWAS)

R_association_testing(analysis_folder=arguments$analysis_folder,
                      phenotype_files=arguments$phenotype_files,
                      covariates=arguments$covariates,
                      GRS_input=arguments$GRS_input,
                      non_GRS_data=arguments$non_GRS_data,
                      group_name_overide=arguments$group_name_overide,
                      PheWAS_manifest_overide=arguments$PheWAS_manifest_overide,
                      analysis_name=arguments$analysis_name,
                      N_cores=arguments$N_cores,
                      phenotype_inclusion_file=arguments$phenotype_inclusion_file,
                      phenotype_exclusion_file=arguments$phenotype_exclusion_file,
                      binary_Case_N=arguments$binary_Case_N,
                      quantitative_Case_N=arguments$quantitative_Case_N)
