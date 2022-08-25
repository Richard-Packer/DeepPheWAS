#! /usr/bin/env Rscript
'This script combines data-field UK Biobank data—collected at the assessment centre—into a single "minimum_tab_data.gz" file that is used throughout the DeepPheWAS pipeline. This is achieved by subsetting the aforementioned data to columns of interest and, where the data are spread across multiple files, merging across all files provided. For UK Biobank data, the column names are referenced by data-field IDs. The desired subset of columns to be extracted is listed in the data_field_ID file found, by default, in data/fields-minimum.txt. This file can be edited to extract a different subset of columns to that provided by default. Files can be inputted in a single folder using the data_folder flag, or specified individually using the data-files flag.

Usage:
    01_minimum_data.R (--data_folder=<FOLDER> | --data_files=<FILES>) (--save_loc=<FILE>) [--r_format --data_field_ID=<FILE> --data_name_pattern=<text> --N_cores=<number> --exclusions=<FILE>]

Options:
    -h --help                     Show this screen.
    -v --version                  Show version.

    Mandated inputs
    --save_loc=<FILE>             Full path to save file location for minimum_tab_data.

    --data_folder=<FOLDER>        Full path of the directory that contains the data files that will be formatted and concatenated.
    Or
    --data_files=<FILES>          Comma separated full file paths of the data files that will be formatted and concatenated.

    Options
    --r_format                    Specific to UK Biobank data. Specify if the input for the data have been downloaded using the R option. If the data
								  have been downloaded using the .csv or .txt options, then no input is required (default).

    --data_field_ID=<FILE>        Full path to the file containing the field_IDs required for Deep-PheWAS if not using default file. Is a plain text
                                  file with no header one field-ID per row. Field-ID is a numeric value, example field-ID 54 is UK biobank assesment
                                  centre.
                                  [default: fields-minimum.txt]

    --data_name_pattern=<text>    Character string for isolating data files if these are in a directory with other files. Defaults to using all files in the directory specified by the data_folder argument.
    --N_cores=<number>            Number of cores requested if parallel computing is desired. Defaults to single core computing.
    --exclusions=<FILE>           Full path to the file containing individuals to be excluded from the analysis. Defaults behaviour is to retain all individuals.

' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.2 01_minimum_data.R')

library(DeepPheWAS)

minimum_data_R(data_folder=arguments$data_folder,
               data_files=arguments$data_files,
               r_format=arguments$r_format,
               data_field_ID=arguments$data_field_ID,
               data_name_pattern=arguments$data_name_pattern,
               N_cores=arguments$N_cores,
               save_loc=arguments$save_loc,
               exclusions=arguments$exclusions)
