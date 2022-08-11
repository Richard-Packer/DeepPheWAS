#! /usr/bin/env Rscript
'
This script combines the results from the split_group split_analysis arguments from 03a_PLINK_assoication_testing.R into a single result file and appends the original results file first removing the empty list and then adding the split group into it.

Usage:
    04_combine_split.R (--split_plink_results_folder=<folder>  --results_RDS_file=<FILE> --group_name=<name>)

Options:
    -h --help  Show this screen.
    -v --version  Show version.

    Mandatory inputs
    --split_plink_results_folder=<folder>       Full path of the folder containing all of the PLINK results from the split_group, split_analysis option in
                                                03a_PLINK_association_testing.R.

    --results_RDS_file=<FILE>                   Full path of the results of the analysis from the 03a_PLINK_association_testing.R script. Is an R object.

    --group_name=<name>                         Name of the group that was split that the results represent.
    --analysis_name=<name>                      Name for the analysis to be used in the naming of the output files, use the same as in 03a_PLINK_association_analysis.

' -> doc


suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.2 04_combine_split.R')

library(devtools)
load_all()

combine_split(split_plink_results_folder=arguments$split_plink_results_folder,
              results_RDS_file=arguments$results_RDS_file,
              group_name=arguments$group_name,
              analysis_name=arguments$analysis_name)
