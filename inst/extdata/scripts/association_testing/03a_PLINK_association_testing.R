#! /usr/bin/env Rscript
'Runs association analysis in using PLINK 2. Takes the variants extracted using extract_snps.R and performs regression analysis on the availble phenotypes per inputted group.

Requires location of a folder where analysis will be hosted, a comma separated list of phenotype files derived from pre_association_preparation.R, a covariate file edited for use in plink (see user guide), and the location of the variants for association prepared with extract_SNPs.R. By default, it will run association tests on all included phenotypes in the phenotype files used and use those phenotype files to assign group names, all analysis is then performed per group. Phenotypes can be specified using one of two arguments phenotype_inclusion_file or phenotype_exclusion_file, group names can be specified using the group_name_overide argument. Results are saved as an R object with option to save tables per-group of the combined raw plink results. With a large number of variants the regression analysis can take a long time, to make for more efficient analysis the phenotypes for any one or more of the analysed groups can be split using the split_group input. When a group is split the phenotype files are split into smaller chunks, the analysis can then be performed in a cluster enviroment.To perfomrm the analysis this same script can be used again but with the split_analysis flag used to amend how results are saved. Results are then combined with 04_combine_split_analysis.R. If splitting analysis tables and figures should not be compiled until 04_combine_split_analysis.R is complete.

Usage:
03a_PLINK_association_testing.R  (--analysis_folder=<FOLDER> --covariate=<file> --phenotype_files=<FILE> --variants_for_association=<FILE> --analysis_name=<name>) [--group_name_overide=<text>] [--split_analysis --PheWAS_manifest_overide=<FILE> --plink_exe=<command> --save_plink_tables --split_group=<name> --N_quant_split=<number> --N_binary_split=<number> --model=<text> --check_existing_results] [--phenotype_inclusion_file=<FILE>  | --phenotype_exclusion_file=<FILE>]

Options:
    -h --help  Show this screen.
    -v --version  Show version.

    Mandatory inputs
    --analysis_folder=<FOLDER>          Full path of the directory in which the association results will be saved.
    --phenotype_files=<FILE>            Comma-separated list containing the full paths of the phenotype data.
    --covariate=<file>                  Full path of a file containing covariate data to be used in model adjustment.
    --variants_for_association=<FILE>   Full path of the genetic data for the SNPs of interest (produced by the extract_SNP.R script).                          
	--analysis_name=<name>              Name for the analysis to be used in the naming of the output files.

    Options
    --group_name_overide=<text>         Comma-separated list containing alternative group names. By default, group names are extracted from 
									    the suffix of the file names provided in the phenotype_files argument. For example, if a file named 
										/home/phenotypes/EUR_phenotypes.csv were provided, the corresponding group name would be "EUR". This 
										argument allows for a different group name to be specified. The order of names provided should match 
										the files specified in the phenotype_files argument.

    --PheWAS_manifest_overide=<FILE>    Full file path of the alternative PheWAS_manifest file.
    --plink_exe=<command>               Command to execute plink2 program [default: plink2]
    --save_plink_tables                 Select if wanting to save the raw plink results per group as a table always saved in
                                        analysis_folder/assoiation_results/group/group_plink_results_raw

    --phenotype_inclusion_file=<FILE>   Full file path to a txt file containing single column containing full PheWAS_ID of phenotypes that will be
                                        included. Cannot be used with phenotype_exclusion_file argument.
    --phenotype_exclusion_file=<FILE>   Full file path to a txt file containing single column containing full PheWAS_ID of phenotypes that will be
                                        excluded. Cannot be used with phenotype_inclusion_file argument.
    --split_group=<name>                Comma separated groups that require splitting for more efficient analysis. Does not run the analysis per group
                                        but splits and saves the phenotypes into smaller files which can then be analysed using cluster based computing.
                                        The number of phenotypes in each split file is dependent on the type (binary or quantitative) and is set using
                                        N_quant_split and N_binary_split. Split files are saved in analysis_folder/group_split with group being the group
                                        name. A file is created in analysis_folder/group_split saved as group_split_guide with group once again being the
                                        group name inputted. This file lists the file names of the splits, which can be used to guide distributed
                                        computing analysis. 
    --N_quant_split=<number>            Number of quantitative phenotypes per split file. [default: 200]
    --N_binary_split=<number>           Number of binary phenotypes per split file. [default: 80]
    --split_analysis                    Specify whether the input being analysed is from the split_group argument.
    --model=<text>                      Genetic model to use for analysis, can be one of genotypic, hethom, dominant, recessive, hetonly. If not one of
                                        these options will be converted to blank.
    --check_existing_results            Input if wanting to re-run part of the analysis whilst checking for existing results in the plink_results folder. Used 
										primarily where a run has aborted part way through.

' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.2 03a_PLINK_association_testing')

library(devtools)
load_all()

plink_association_testing(analysis_folder=arguments$analysis_folder,
                          phenotype_files=arguments$phenotype_files,
                          covariate=arguments$covariate,
                          variants_for_association=arguments$variants_for_association,
                          analysis_name=arguments$analysis_name,
                          group_name_overide=arguments$group_name_overide,
                          PheWAS_manifest_overide=arguments$PheWAS_manifest_overide,
                          plink_exe=arguments$plink_exe,
                          save_plink_tables=arguments$save_plink_tables,
                          phenotype_inclusion_file=arguments$phenotype_inclusion_file,
                          phenotype_exclusion_file=arguments$phenotype_exclusion_file,
                          split_group=arguments$split_group,
                          N_quant_split=arguments$N_quant_split,
                          N_binary_split=arguments$N_binary_split,
                          split_analysis=arguments$split_analysis,
                          model=arguments$model,
                          check_existing_results=arguments$check_existing_results)
