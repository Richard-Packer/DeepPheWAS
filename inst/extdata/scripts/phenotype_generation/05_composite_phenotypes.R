#!/usr/bin/env Rscript
'This script creates composite phenotypes and composite concepts using a central map file. Also creates population_control ID lists that are used when creating the composite phenotypes. The phenotype creating script iterates over itself (default five times) as some composite phenotypes use other composite phenotypes in their definition.

Usage:
    05_composite_phenotypes.R (--phenotype_save_file=<FILE>) (--phenotype_folder=<FOLDER> | --phenotype_files=<FILES>)[--control_populations=<FILE> --N_iterations=<number> --update_list=<FILE> --composite_phenotype_map_overide=<FILE> --control_pop_save_file=<FILE>]

Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.

    Mandatory
    --phenotype_save_file=<FILE>              Full path for the save file for the generated composite phenotypes RDS.
    --phenotype_folder=<FOLDER>               Full path of the folder containing the phenotype data created in previous steps.
    --phenotype_files=<FILES>                 Comma separated full file paths of phenotype data created in previous steps.

    Options
    --control_pop_save_file=<FILE>            Full path for the save file for the generated composite control populations RDS.
    --composite_phenotype_map_overide=<FILE>  Full path of the composite_phenotyope_map file. Provided with the package.

    --N_iterations=<number>                   Number of iterations curated phenotype script is to run through. [default: 5]

    --control_populations=<FILE>              Full path of the "control_populations" file containing columns of IDs, column names are used to create lists of IDs.
											  Saved as a list called control_populations.RDS in the folder inputted in the phenotype_folder flag. The lists are used to
											  define some of the composite phenotype control populations as directed by the composite_phenotype_map file. The unedited
									          version of the composite_phenotype_map file uses two populations, all_pop and primary_care_pop, which represent all IDs
											  in the sample and all IDs with available primary care data. To create composite phenotypes, the names of these list must match
											  the composite_phenotype_map file. The 02_data_preparation.R script creates a control_populations file that is used by default.

    --update_list=<FILE>                      Option to run this script as an update to an existing composite_phenotype list object. If specified, this script will use
											  phenotype_save_file to load the existing list and then save over that file upon completion.
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.2 05_composite_phenotypes.R')

library(DeepPheWAS)

composite_phenotyping(composite_phenotype_map_overide=arguments$composite_phenotype_map_overide,
                      phenotype_folder=arguments$phenotype_folder,
                      phenotype_files=arguments$phenotype_files,
                      N_iterations=arguments$N_iterations,
                      phenotype_save_file=arguments$phenotype_save_file,
                      control_populations=arguments$control_populations,
                      update_list=arguments$update_list,
                      control_pop_save_file=arguments$control_pop_save_file)
