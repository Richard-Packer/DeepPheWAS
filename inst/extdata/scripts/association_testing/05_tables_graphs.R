#!/usr/bin/env Rscript
'Creates result tables and graphs for the association analyses from either the GRS or plink methods. The two usages below represent inputting data from 03a_PLINK_association_testing.R and 03b_R_association_testing.R.
To produce graphs input the per_group_name_graph, per_snp_graph or R_association_graph depending on result source and transformation of the data required. There are many options to select related to filtering for MAC in plink results and general appearance of the graphs. By using the provided filter inputs it is possible to edit any individual graphs using the array of options.

Usage: 05_tables_graphs.R (--results_file=<FILE> --analysis_name=<name> --plink_results --SNP_list=<FILE> --save_folder=<FOLDER>) [ --group_filter=<text> --PheWAS_ID_filter=<FILE> --PheWAS_manifest_overide=<FILE> --max_pheno=<number> --sig_FDR=<number> --no_save_table_all_results --no_graph_all --no_graph_sig --max_FDR_graph=<number> --SNP_filter=<FILE> --group_name_filter=<FILE> --save_raw_plink --MAC=<number> --MAC_case=<number> --MAC_control=<number> --per_group_name_graph --per_snp_graph --save_table_per_group_name --save_table_per_snp --sex_split --PheWAS_ID_label_filter=<FILE> --max_overlap_labels=<number> --graph_file_save=<name> --label_text_size=<number> --order_groups_alphabetically --order_phenotypes_alphabetically --save_all_graphs]

       05_tables_graphs.R (--results_file=<FILE> --analysis_name=<name> --R_association_results --save_folder=<FOLDER>) [--group_filter=<text> --PheWAS_ID_filter=<FILE> --PheWAS_manifest_overide=<FILE> --max_pheno=<number> --sig_FDR=<number> --no_save_table_all_results --no_graph_all --no_graph_sig --max_FDR_graph=<number> --R_association_graph --sex_split --PheWAS_ID_label_filter=<FILE> --max_overlap_labels=<number> --graph_file_save=<name> --label_text_size=<number> --order_groups_alphabetically --order_phenotypes_alphabetically --save_all_graphs]

Options:
    -h --help  Show this screen.
    -v --version  Show version.

    Mandatory with any input
    --results_file=<FILE>             Full file path of the results file RDS R list object.
    --analysis_name=<name>            Name for the analysis.
	  --save_folder=<FOLDER>            Full file path of the folder to which the output will be saved.

    Plink result input
    --plink_results                   Select if results are from 03a_PLINK_association_testing.R
    --SNP_list=<FILE>                 Full path of the SNP_list file used in 02_extracting_snps.R.

    R association testing input
    --R_association_results           Select if results are from 03b_R_association_testing.R

    Options for both inputs

    --group_filter=<text>               Comma-separated text input, used to filter the group to which the table and graph functions are applied. Inputted groups
									                      are the ones that are retained for analysis, group here refers to the grouping variable used to subset the analysis classically ancestry.
    --PheWAS_ID_filter=<FILE>           Full path of a file describing the subset of PheWAS IDs that will be included in the output.
    --PheWAS_manifest_overide=<FILE>    Full file path of the alternative PheWAS_manifest file.
    --max_pheno=<number>                Manual override for inputting maximum phenotypes analysed. Used for calculating FDR. The default used the largest number
                                        of associations in the results_file of any grouping.
    --sig_FDR=<number>                  Value of FDR for which associations are reported as significant. Will alter the line of significance in both graph types
                                        and the reported phenotypes in the significant_pheno graphs. [default: 0.01]
    --no_save_table_all_results         Specify whether a table containing all the results should be saved (FALSE) or not (TRUE). Default is FALSE.
    --no_graph_all                      Specify whether a figure of all the results should be produced (FALSE) or not (TRUE). Default is FALSE.
    --no_graph_sig                      Specify whether a figure of the significant results should be produced (FALSE) or not (TRUE). Default is FALSE.
    --max_FDR_graph=<number>            Maximum value for the y-axis to be used when the FDR for a particular association is exceptionally small. Default is 300. [default: 300]
    --sex_split                         Specify whether the figure should be stratified by sex (TRUE) or not (FALSE). If TRUE, this will produce separate figures for
                                        males, females and combined. Default is FALSE. Does not create split tables.
    --PheWAS_ID_label_filter=<FILE>     Full file path to plain text file containing single headed column of PheWAS_IDs. Only these PheWAS_IDs will be labelled within any graphical
                                        output designed to be used primarily when trying to edit a single graph as the filter will apply to all graphs being created.
    --max_overlap_labels=<number>       Number, represents maximum overlaps for labelling of phenotypes in the all_pheno graph, lowering the number has the effect of reducing the total
                                        number of phenotypes labelled. [default: 20]
    --graph_file_save=<name>            Allows user to specify the file format of the graphs, for example pdf or png. [default: png]
    --label_text_size=<number>          Number, represents the text size of the labelled phenotypes in all_pheno graph. [default: 2]
    --order_groups_alphabetically       Specify whether the groups in the graphs should be ordered by lowest FDR (FALSE) or alphabetically (TRUE).
    --order_phenotypes_alphabetically   Specify whether the phenotypes within each group are ordered by the lowest FDR (FALSE) or alphabetically (TRUE).
    --save_all_graphs                   Specify whether to always save every graph with or without a significant result (TRUE) or to only save when at least one association is
                                        significant (FALSE).

    Options for plink results
    --SNP_filter=<FILE>                 Full path of a file containing SNP IDs, used to filter the output to a given subset of SNPs. Use if wanting to apply the table and/or graphing functions
                                        a subset of results. Only works for plink_results.
    --group_name_filter=<FILE>          Full path of the file containing a list of group_names that correspond with those defined in SNP_list. This can be used if the table and/or
									                      graphing functions should be applied to a subset of SNPs, such as only those in a given credible set. This option can only be used with
									                      association results generated by 03a_PLINK_association_testing.R.
    --save_raw_plink                    Specify whether the unfiltered association results from the PLINK analysis should be saved (TRUE) or not (FALSE). Default is FALSE.
    --MAC=<number>                      MAC (minor allele count) filter applied to all associations, only applicable in results from03a_PLINK_association_testing.R
                                        [default: 20]
    --MAC_case=<number>                 Minor allele count among cases for filtering the association results. only applicable in results from 03a_PLINK_association_testing.R and only in
                                        binary phenotypes. Default is 5. [default: 5]
    --MAC_control=<number>              MAC (minor allele count) filter applied only to controls, only applicable in results from PheWAS_association_PLINK.R
                                        and only in binary phenotypes. [default: 10]
    --per_group_name_graph              Select if wanting to produce graphs per-group_name. This is used when looking to report the most significant finding
                                        across several SNPs for a single construct, potentially and gene or a sentinal SNP with a credible set. It is an
                                        column in the SNP_list file. Default is FALSE.
    --per_snp_graph                     Specify whether a figure is produced for every SNP provided (TRUE) or not (FALSE). Default is FALSE.
    --save_table_per_group_name         Select if wanting to save a table of results per_group_name. Will extract the lowest FDR value for a
                                        PheWAS_ID-group_name combination and report that. Only works for plink_results.
    --save_table_per_snp                Specify whether a results table should be generated for every SNP provided (TRUE) or not (FALSE). Default is FALSE. Will be saved in a created folder named
                                        /analysis_name_group_per_SNP_tables. Example if the group was groupA and analysis_name top_SNPs the folder would be
                                        /top_SNPs_groupA_per_SNP_tables.

    Options R association results
    --R_association_graph             Specify whether a figure of the results from the association analysis from 03b_R_association_testing.R should be produced (TRUE) or not (FALSE). Default is FALSE.
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc)

library(DeepPheWAS)

graphs_tables_DeepPheWAS(results_file=arguments$results_file,
                         analysis_name=arguments$analysis_name,
                         plink_results=arguments$plink_results,
                         SNP_list=arguments$SNP_list,
                         R_association_results=arguments$R_association_results,
                         save_folder=arguments$save_folder,
                         group_filter=arguments$group_filter,
                         PheWAS_ID_filter=arguments$PheWAS_ID_filter,
                         PheWAS_manifest_overide=arguments$PheWAS_manifest_overide,
                         max_pheno=arguments$max_pheno,
                         sig_FDR=arguments$sig_FDR,
                         no_save_table_all_results=arguments$no_save_table_all_results,
                         no_graph_all=arguments$no_graph_all,
                         no_graph_sig=arguments$no_graph_sig,
                         max_FDR_graph=arguments$max_FDR_graph,
                         sex_split=arguments$sex_split,
                         SNP_filter=arguments$SNP_filter,
                         group_name_filter=arguments$group_name_filter,
                         save_raw_plink=arguments$save_raw_plink,
                         MAC=arguments$MAC,
                         MAC_case=arguments$MAC_case,
                         MAC_control=arguments$MAC_control,
                         per_group_name_graph=arguments$per_group_name_graph,
                         per_snp_graph=arguments$per_snp_graph,
                         save_table_per_group_name=arguments$save_table_per_group_name,
                         save_table_per_snp=arguments$save_table_per_snp,
                         R_association_graph=arguments$R_association_graph,
                         PheWAS_ID_label_filter=arguments$PheWAS_ID_label_filter,
                         max_overlap_labels=arguments$max_overlap_labels,
                         graph_file_save=arguments$graph_file_save,
                         label_text_size=arguments$label_text_size,
                         order_groups_alphabetically=arguments$order_groups_alphabetically,
                         order_phenotypes_alphabetically=arguments$order_phenotypes_alphabetically,
                         save_all_graphs=arguments$save_all_graphs
                         )
