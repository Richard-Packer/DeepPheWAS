#! /usr/bin/env Rscript
'This script creates quantitative phenotypes from UK Biobank primary care data. The script uses the PheWAS_manifest.csv file (provided with package) as a guide to creating phenotypes.

Usage:
    03d_primary_care_quantititive_phenotypes.R (--GPC=<FILE> --DOB=<FILE> --phenotype_save_file=<FILE>) [--N_cores=<number> --PheWAS_manifest_overide=<FILE> --code_list_overide=<FOLDER>]

Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.

    Mandated inputs
    --GPC=<FILE>                              Full path of the primary care clinical data file from UK Biobank.
    --DOB=<FILE>                              Full path of a file describing participant date of birth. This does not need to be the exact date. By default for UK Biobank data,
											  the script creates the date of birth from month and year of birth (data-fields 52 and 34).
    --phenotype_save_file=<FILE>              Full path of the save file for the generated concepts RDS to be used for phenotype creation.

    Options
    --N_cores=<number>                        Number of cores requested if parallel computing is desired. Defaults to single core computing.
    --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.
    --code_list_overide=<FOLDER>              Full file path of the folder containing alternative primary care code lists. [default: default]

' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.2 03d_primary_care_quantitiative_phenotypes.R')

library(DeepPheWAS)

primarycare_quantitative_phenotypes(GPC=arguments$GPC,
                                    DOB=arguments$DOB,
                                    phenotype_save_file=arguments$phenotype_save_file,
                                    N_cores=arguments$N_cores,
                                    PheWAS_manifest_overide=arguments$PheWAS_manifest_overide,
                                    code_list_overide=arguments$code_list_overide)
