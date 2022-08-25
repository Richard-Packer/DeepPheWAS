#! /usr/bin/env Rscript
'This script uses code lists provided to define concepts. Concepts can then be combined to form cases and controls for "composite-phenotypes". A concept is a single homogenous group of codes that all describe a single disease, symptom, or related group of medicines.
The main output that is used for forming phenotypes is saved as an R object and is a record of the ID, number of codes and the date of the earliest code in a record. The other files are used for creating summary documents of the concepts.
UK Biobank prescription data is unusual in that it does not contain consistent codes (namely BNF codes) for most of the data. As such, much of the searching is done using key search terms. A small subsection requires using Read version 2 drug codes. Because of this major difference with the other types of data, prescription concepts have their own function and would need to be used with caution in a non-UK Biobank setting.
Users can choose to define only concepts that use clinical data or prescription data or both. At least one must be inputted, however.

Usage:
    03c_creating_concepts.R (--GPP=<FILE> --concept_save_file=<FILE> --all_dates_save_file=<FILE>) [--health_data=<FILE>  --PheWAS_manifest_overide=<FILE> --code_list_folder_override=<FOLDER>]
    03c_creating_concepts.R (--health_data=<FILE> --concept_save_file=<FILE> --all_dates_save_file=<FILE>) [--GPP=<FILE> --PheWAS_manifest_overide=<FILE> --code_list_folder_override=<FOLDER>]

Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.

    Mandatory inputs
    --concept_save_file=<FILE>                Full path of the save file for the generated concepts RDS to be used for phenotype creation.
    --all_dates_save_file=<FILE>              Full file path for the save file for the generated all_dates.RDS used for per-event combinations of concepts.

    Options
    --health_data=<FILE>                      Full path of the health_data file for UK Biobank. Specifying this argument will generate clinical concepts.
    --GPP=<FILE>                              Full path of the primary care prescription data from UK Biobank. Specifying this argument will generate prescription concepts.

    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
    --code_list_folder_override=<FOLDER>      Full file for the folder containing code lists only use if not using default stored in R package.

' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.2 03c_creating_concepts.R')

library(DeepPheWAS)

creating_concepts(concept_save_file=arguments$concept_save_file,
                  all_dates_save_file=arguments$all_dates_save_file,
                  health_data=arguments$health_data,
                  GPP=arguments$GPP,
                  PheWAS_manifest_overide=arguments$PheWAS_manifest_overide,
                  code_list_folder_override=arguments$code_list_folder_override)
