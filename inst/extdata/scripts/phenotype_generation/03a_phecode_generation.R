#!/usr/bin/env Rscript
'This script maps healthcare data from the ICD-9 and ICD-10 coding frameworks to Phecodes (https://phewascatalog.org/phecodes) and generates output that describes the resulting phenotypes in all individuals. In addition, output can be generated that describes individuals who should be excluded when selecting controls (binary phenotypes). For example, if one were interested in defining an asthma phenotype but wanted any controls to not have chronic obstructive pulmonary disease codes, that can be done using this script.

Usage:
    03a_phecode_generation.R (--health_data=<FILE> --sex_info=<FILE>) (--phecode_save_file=<FILE> | --no_phecodes) (--no_range_ID | --range_ID_save_file=<FILE>) [--control_exclusions=<FILE> --N_cores=<number> --ICD10=<text_comma> --ICD9=<text_comma> --PheWAS_manifest_overide=<FILE>]

Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.

    Mandated inputs
    --health_data=<FILE>                      Full path of the file produced by the previous step (phenotype_data.R).

    --sex_info=<FILE>                         Full file path of the combined_sex file a file containing participant ID and sex information
                                              (0=female, 1=male).

    Use at least one of
    --no_phecodes                             Choose whether to save the Phecode-generated phenotype data (FALSE) or not (TRUE). Defaults to FALSE.
    OR
    --phecode_save_file=<FILE>                Full path of the save file for the Phecode phenotypes R data object.

    Use at least one of      z
    --range_ID_save_file=<FILE>               Full path of the folder used to store the control exclusions R data object.

    Options
    --control_exclusions=<FILE>               Full file path of the optional control exclusions file.
    --N_cores=<number>                        Number of cores requested if parallel computing is desired. Defaults to single core computing.
    --no_range_ID                             Choose whether to save lists of participant identifiers to be excluded from control definitions (FALSE) or not (TRUE). Defaults to FALSE.

    --ICD10=<text_comma>                      Comma-separated string representing the column names used in health_data for ICD-10 values. If there are no ICD-10 values,
											                        use NA or any text not used as a source in health data. Defaults to "ICD10".

    --ICD9=<text_comma>                       Comma-separated string representing the column names used in health_data for ICD-9 values. If there are
											                        no ICD-9 values, use NA or any text not used as a source in health data. Defaults to "ICD9".
	  --PheWAS_manifest_overide=<FILE>          Full file path of the alternative PheWAS_manifest file.


' -> doc

suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'v0.1 pheocode_generation.R')

library(DeepPheWAS)

generating_phecodes(health_data=arguments$health_data,
                    sex_info=arguments$sex_info,
                    phecode_save_file=arguments$phecode_save_file,
                    range_ID_save_file=arguments$range_ID_save_file,
                    control_exclusions=arguments$control_exclusions,
                    N_cores=arguments$N_cores,
                    no_range_ID=arguments$no_range_ID,
                    no_phecodes=arguments$no_phecodes,
                    ICD10=arguments$ICD10,
                    ICD9=arguments$ICD10,
                    PheWAS_manifest_overide=arguments$PheWAS_manifest_overide)
