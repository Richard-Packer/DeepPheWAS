#' Produces graphs for R-association data.
#'
#' @param x Data table for graphs
#' @param save_root full file location for save file minus specific file name
#' @param analysis_name name of teh analysis
#' @param graph_choice numeric value indicating if all or only some graphs are chosen
#' @param FDR_figure value of significance for FDR.
#' @param max_FDR_graph value for lowest fdr value, converts 0 into this number.
#' @param PheWAS_label_filter list of PheWAS_IDs to label within the all_pheno graph.
#' @param max_overlap maximum number of overlaps for labelled phenotypes in the all_pheno graphs
#' @param graph_type save format of the graphs any input readable from ggsave is accepted
#' @param label_size size of the text for labelled phenotypes.
#' @param order_groups_alphabetically T or F to order groups on graph alphabetically rather than default lowest FDR.
#' @param order_phenotypes_alphabetically T or F to order phenotypes on graphs alphabetically rather than default lowest FDR.
#' @param save_all_graphs T or F to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return saved graphs per group
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
R_association_function <- function(x,save_root,analysis_name,graph_choice,FDR_figure,max_FDR_graph,PheWAS_label_filter,max_overlap,graph_type,label_size,order_groups_alphabetically,
                                   order_phenotypes_alphabetically,
                                   save_all_graphs){
  R_association_graph <- x %>%
    dplyr::group_split(.data$name)
  sex_label <- unique(x$sex_pheno_identifier)
  graph_save_location <- save_root

  mapply(making_graphs,R_association_graph,MoreArgs = list(b="R",c=graph_choice,d=sex_label,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,graph_save_location=graph_save_location,PheWAS_label_filter=PheWAS_label_filter,max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                                           order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                           save_all_graphs=save_all_graphs))
}

#' Produces per-group graphs.
#'
#' @param x Data table for graphs
#' @param y Group name
#' @param save_root full file location for save file minus specific file name
#' @param analysis_name name of teh analysis
#' @param graph_choice numeric value indicating if all or only some graphs are chosen
#' @param FDR_figure value of significance for FDR.
#' @param max_FDR_graph value for lowest fdr value, converts 0 into this number.
#' @param PheWAS_label_filter list of PheWAS_IDs to label within the all_pheno graph.
#' @param max_overlap maximum number of overlaps for labelled phenotypes in the all_pheno graphs
#' @param graph_type save format of the graphs any input readable from ggsave is accepted
#' @param label_size size of the text for labelled phenotypes.
#' @param order_groups_alphabetically T or F to order groups on graph alphabetically rather than default lowest FDR.
#' @param order_phenotypes_alphabetically T or F to order phenotypes on graphs alphabetically rather than default lowest FDR.
#' @param save_all_graphs T or F to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return saved graphs per group
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
per_group_function <- function(x,y,save_root,analysis_name,graph_choice,FDR_figure,max_FDR_graph,PheWAS_label_filter,max_overlap,graph_type,label_size,order_groups_alphabetically,
                               order_phenotypes_alphabetically,
                               save_all_graphs) {
  per_group_name <- x %>%
    dplyr::group_by(.data$collective_name,.data$PheWAS_ID) %>%
    dplyr::slice_min(.data$FDR) %>%
    dplyr::ungroup() %>%
    dplyr::group_split(.data$collective_name)

  per_group_name_folder <- paste0(save_root,"/",analysis_name,"_",y,"_results_per_collective_name/")
  dir.create(per_group_name_folder)

  graph_save_location <- per_group_name_folder

  sex_label <- unique(x$sex_pheno_identifier)

  mapply(making_graphs,per_group_name,MoreArgs = list(b="per_group_name",c=graph_choice,d=sex_label,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,graph_save_location=graph_save_location,PheWAS_label_filter=PheWAS_label_filter,max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                                      order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                      save_all_graphs=save_all_graphs))
}
#' Produces per-SNP graphs.
#'
#' @param a Data table for graphs
#' @param y Group name
#' @param save_root full file location for save file minus specific file name
#' @param analysis_name name of the analysis
#' @param graph_choice numeric value indicating if all or only some graphs are chosen
#' @param FDR_figure value of significance for FDR.
#' @param max_FDR_graph value for lowest fdr value, converts 0 into this number.
#' @param PheWAS_label_filter list of PheWAS_IDs to label within the all_pheno graph.
#' @param max_overlap maximum number of overlaps for labelled phenotypes in the all_pheno graphs
#' @param graph_type save format of the graphs any input readable from ggsave is accepted
#' @param label_size size of the text for labelled phenotypes.
#' @param order_groups_alphabetically T or F to order groups on graph alphabetically rather than default lowest FDR.
#' @param order_phenotypes_alphabetically T or F to order phenotypes on graphs alphabetically rather than default lowest FDR.
#' @param save_all_graphs T or F to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return saved graphs per SNP.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
per_snp_function <- function(a,y,save_root,analysis_name,graph_choice,FDR_figure,max_FDR_graph,PheWAS_label_filter,max_overlap,graph_type,label_size,order_groups_alphabetically,
                             order_phenotypes_alphabetically,
                             save_all_graphs){
  per_SNP <- a %>%
    dplyr::group_split(.data$ID)
  sex_label <- unique(a$sex_pheno_identifier)
  per_SNP_folder <- paste0(save_root,"/",analysis_name,"_",y,"_results_per_SNP/")
  dir.create(per_SNP_folder)
  graph_save_location <- per_SNP_folder

  mapply(making_graphs,per_SNP,MoreArgs = list(b="per_snp",c=graph_choice,d=sex_label,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,graph_save_location=graph_save_location,PheWAS_label_filter=PheWAS_label_filter,max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                               order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                               save_all_graphs=save_all_graphs))
}
#' Converts P values into false discovery rate.
#'
#' @param x Group name
#' @param max_pheno_tests the maximum number of tests that are being run across grouping variable. Must be at => the number of associations.
#' @return data frame with new fdr column.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data

fdr_calc <- function(x,max_pheno_tests) {
  fdr <- x %>%
    dplyr::distinct(.data$PheWAS_ID, .keep_all = T) %>%
    dplyr::mutate(FDR=stats::p.adjust(p = .data$P, method = "fdr", n = max_pheno_tests)) %>%
    dplyr::arrange(.data$FDR)

  return(fdr)
}
#' Makes graphs from converted RDS files.
#'
#' @param a type of analysis data
#' @param b Type of analysis
#' @param c Graph choice
#' @param d sex phenotype labels
#' @param FDR_figure value for significance threshold for fdr.
#' @param max_FDR_graph value for lowest fdr value, converts 0 into this number.
#' @param graph_save_location full file location to save graphs
#' @param PheWAS_label_filter list of PheWAS_IDs to label within the all_pheno graph
#' @param max_overlap maximum number of overlaps for labelled phenotypes in the all_pheno graphs
#' @param graph_type save format of the graphs any input readable from ggsave is accepted
#' @param label_size size of the text for labelled phenotypes.
#' @param order_groups_alphabetically T or F to order groups on graph alphabetically rather than default lowest FDR.
#' @param order_phenotypes_alphabetically T or F to order phenotypes on graphs alphabetically rather than default lowest FDR.
#' @param save_all_graphs T or F to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return Creates and then saves graphs from results files.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
making_graphs <- function(a,b,c,d,FDR_figure,max_FDR_graph,graph_save_location,PheWAS_label_filter,max_overlap,graph_type,label_size,order_groups_alphabetically,
                          order_phenotypes_alphabetically,
                          save_all_graphs) {
  . <- NULL
  data <- a
  data_test <- data %>%
    dplyr::filter(.data$FDR<=FDR_figure)
  if(isFALSE(save_all_graphs) & nrow(data_test)<1){
    return()
  } else {
    data <- data %>%
      dplyr::mutate(FDR=ifelse(.data$FDR==0,max_FDR_graph,.data$FDR))

    if(order_groups_alphabetically){
      order <- data  %>%
        dplyr::group_by(.data$phenotype_group) %>%
        dplyr::summarise(min_FDR=min(.data$FDR,na.rm = T)) %>%
        dplyr::arrange(.data$phenotype_group) %>%
        dplyr::mutate(group_number =seq(1:(nrow(.)))) %>%
        dplyr::select(.data$phenotype_group,.data$group_number)
    } else {
      order <- data  %>%
        dplyr::group_by(.data$phenotype_group) %>%
        dplyr::summarise(min_FDR=min(.data$FDR,na.rm = T)) %>%
        dplyr::arrange(.data$min_FDR) %>%
        dplyr::mutate(group_number =seq(1:(nrow(.)))) %>%
        dplyr::select(.data$phenotype_group,.data$group_number)
    }
    if(order_phenotypes_alphabetically){
      if(is.null(PheWAS_label_filter)){
        all_pheno_data <- data %>%
          dplyr::left_join(order,by="phenotype_group") %>%
          tidyr::drop_na(.data$FDR) %>%
          dplyr::arrange(.data$description) %>%
          dplyr::arrange(.data$group_number) %>%
          dplyr::mutate(annotate=ifelse(.data$FDR<=FDR_figure,T,F),
                        short_desc_T=paste0(.data$short_desc," ",.data$group),
                        short_desc=factor(.data$short_desc, levels = unique(.$short_desc)),
                        seq=seq(nrow(.)),
                        y_axis_info=-log10(.data$FDR)) %>%
          dplyr::mutate(short_desc_T=factor(.data$short_desc_T, levels = unique(.$short_desc_T)),
                        seq=seq+150*(.data$group_number-1))
      } else {
        all_pheno_data <- data %>%
          dplyr::left_join(order,by="phenotype_group") %>%
          tidyr::drop_na(.data$FDR) %>%
          dplyr::arrange(.data$description) %>%
          dplyr::arrange(.data$group_number) %>%
          dplyr::mutate(annotate=ifelse(.data$PheWAS_ID %in% PheWAS_label_filter,T,F),
                        short_desc_T=paste0(.data$short_desc," ",.data$group),
                        short_desc=factor(.data$short_desc, levels = unique(.$short_desc)),
                        seq=seq(nrow(.)),
                        y_axis_info=-log10(.data$FDR)) %>%
          dplyr::mutate(short_desc_T=factor(.data$short_desc_T, levels = unique(.$short_desc_T)),
                        seq=seq+150*(.data$group_number-1))
      }

    } else {

      if(is.null(PheWAS_label_filter)){
        all_pheno_data <- data %>%
          dplyr::left_join(order,by="phenotype_group") %>%
          tidyr::drop_na(.data$FDR) %>%
          dplyr::arrange(.data$FDR) %>%
          dplyr::arrange(.data$group_number) %>%
          dplyr::mutate(annotate=ifelse(.data$FDR<=FDR_figure,T,F),
                        short_desc_T=paste0(.data$short_desc," ",.data$group),
                        short_desc=factor(.data$short_desc, levels = unique(.$short_desc)),
                        seq=seq(nrow(.)),
                        y_axis_info=-log10(.data$FDR)) %>%
          dplyr::mutate(short_desc_T=factor(.data$short_desc_T, levels = unique(.$short_desc_T)),
                        seq=seq+150*(.data$group_number-1))
      } else {
        all_pheno_data <- data %>%
          dplyr::left_join(order,by="phenotype_group") %>%
          tidyr::drop_na(.data$FDR) %>%
          dplyr::arrange(.data$FDR) %>%
          dplyr::arrange(.data$group_number) %>%
          dplyr::mutate(annotate=ifelse(.data$PheWAS_ID %in% PheWAS_label_filter,T,F),
                        short_desc_T=paste0(.data$short_desc," ",.data$group),
                        short_desc=factor(.data$short_desc, levels = unique(.$short_desc)),
                        seq=seq(nrow(.)),
                        y_axis_info=-log10(.data$FDR)) %>%
          dplyr::mutate(short_desc_T=factor(.data$short_desc_T, levels = unique(.$short_desc_T)),
                        seq=seq+150*(.data$group_number-1))
      }
    }
    if(all(!is.na(match(c("Beta","OR"),colnames(all_pheno_data))))){
      all_pheno_data <- all_pheno_data %>%
        dplyr::mutate(size=dplyr::case_when(is.na(.data$OR) ~ sqrt(.data$Beta^2),
                                            is.na(.data$Beta) ~ sqrt((log(.data$OR))^2)))
    } else if(all(!is.na(match("Beta",colnames(all_pheno_data))) && is.na(match("OR",colnames(all_pheno_data))))) {
      all_pheno_data <- all_pheno_data %>%
        dplyr::mutate(size=sqrt(.data$Beta^2))
    } else if(all(!is.na(match("OR",colnames(all_pheno_data))) && is.na(match("Beta",colnames(all_pheno_data))))) {
      all_pheno_data <- all_pheno_data %>%
        dplyr::mutate(size=sqrt((log(.data$OR))^2))
    }


    labels= dplyr::summarize(dplyr::group_by(all_pheno_data, .data$group_number), tick=mean(unique(seq)),label=as.character(.data$phenotype_group[1]))
    labels=labels[order(labels$tick),]

    round_any = function(x, accuracy, f=ceiling){f(x/ accuracy) * accuracy}

    shapes = c("positive" = 24, "negative" = 25)
    max.x<-max(all_pheno_data$seq)
    max.y <- round_any(max(all_pheno_data$y_axis_info),5)

    if(b=="per_group_name") {
      name_file <- paste0(unique(all_pheno_data$collective_name),d)
    } else if(b=="per_snp") {
      name_file <- paste0(unique(all_pheno_data$graph_save_name),"_",unique(all_pheno_data$group),d)
    } else if(b=="R") {
      name_file <- paste0(unique(all_pheno_data$name_group),d)
    }


    # making the graphs
    all_pheno_graph <- ggplot2::ggplot(all_pheno_data, ggplot2::aes(x = seq, y = -log10(.data$FDR), color = .data$phenotype_group, fill = .data$phenotype_group)) +
      ggplot2::geom_point(ggplot2::aes(shape=.data$effect_direction)) +
      ggplot2::scale_shape_manual(values = shapes) +
      viridis::scale_colour_viridis(option = "H",discrete=T) +
      viridis::scale_fill_viridis(option= "H",discrete=T) +
      ggplot2::scale_y_continuous(limits=c(0,max.y)) +
      ggplot2::scale_x_continuous(name="Phenotype groups", limits=c(1,max.x), breaks=labels$tick, labels=labels$label, expand=c(.01,0)) +
      ggplot2::geom_hline(yintercept=-log10(FDR_figure),colour="red", alpha=I(1/3),size=1) +
      ggrepel::geom_text_repel(ggplot2::aes(label=.data$short_desc),colour="black",data=all_pheno_data[all_pheno_data$annotate,],
                               size=label_size,angle=0,max.overlaps = max_overlap) +
      ggplot2::guides(color="none",fill="none",shape=ggplot2::guide_legend("Direction of Effect")) +
      ggplot2::scale_size_area(name="Effect size",breaks=scales::breaks_extended(3),max_size = 4) +
      ggplot2::theme(panel.background=ggplot2::element_blank(),
                     panel.grid.minor=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=9, colour="black", hjust=1, vjust=.5),
                     axis.text.x=ggplot2::element_text(size=7, colour="black", angle=-45, hjust=0, vjust=0),
                     axis.line =ggplot2::element_line(colour="black"),
                     axis.ticks=ggplot2::element_line(colour="black"),
                     legend.position = c(0.8,0.8),
                     legend.key=ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(1,60,1,1))

    sig_pheno_data <- all_pheno_data %>%
      dplyr::filter(.data$FDR<0.01)

    sig_results_all <- ggplot2::ggplot(sig_pheno_data, ggplot2::aes(x = .data$short_desc, y = -log10(.data$FDR), color = .data$phenotype_group, fill = .data$phenotype_group)) +
      ggplot2::geom_point(ggplot2::aes(shape=.data$effect_direction)) +
      ggplot2::scale_shape_manual(values = shapes) +
      ggplot2::scale_x_discrete(name="Phenotypes") +
      viridis::scale_colour_viridis(option = "H",discrete=T) +
      viridis::scale_fill_viridis(option= "H",discrete=T) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=-log10(FDR_figure)),colour="red", alpha=I(1/3),size=1) +
      ggplot2::scale_y_continuous(limits=c(0, max(-log10(sig_pheno_data$FDR)))) +
      ggplot2::guides(fill=ggplot2::guide_legend(title="Phenotypic Category", ncol=1, override.aes=list(shape=24,size=1)),
                      col=ggplot2::guide_legend(title="Phenotypic Category", ncol=1),
                      shape=ggplot2::guide_legend(title="Direction of Effect", ncol=1,override.aes=list(size=1))) +
      ggplot2::scale_size_area(name="Effect size",breaks=scales::breaks_extended(3),max_size = 4) +
      ggplot2::theme(panel.background=ggplot2::element_blank(),
                     panel.grid.minor=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_text(size=9, colour="black", hjust=1, vjust=.5),
                     axis.text.x=ggplot2::element_text(size=6, colour="black", angle=90, hjust=1, vjust=0),
                     axis.line =ggplot2::element_line(colour="black"),
                     axis.ticks=ggplot2::element_line(colour="black"),
                     legend.key=ggplot2::element_blank(),
                     legend.key.height=ggplot2::unit(0.5,"line"))
    if(nrow(data_test)<1){
      ggplot2::ggsave(filename = paste0(graph_save_location,"/",name_file,"_all_pheno.",graph_type),
                      plot = all_pheno_graph, device = graph_type, dpi = 300, width = 9,height = 6,units = "in")
    } else if(c=="both") {
      ggplot2::ggsave(filename = paste0(graph_save_location,"/",name_file,"_sig_pheno.",graph_type),
                      plot = sig_results_all, device = graph_type, dpi = 300, width = 9,height = 6,units = "in")
      ggplot2::ggsave(filename = paste0(graph_save_location,"/",name_file,"_all_pheno.",graph_type),
                      plot = all_pheno_graph, device = graph_type, dpi = 300, width = 9,height = 6,units = "in")

    } else if (c=="sig_only") {
      ggplot2::ggsave(filename = paste0(graph_save_location,"/",name_file,"_sig_pheno.",graph_type),
                      plot = sig_results_all, device = graph_type, dpi = 300, width = 9,height = 6,units = "in")

    } else if(c=="all_pheno") {
      ggplot2::ggsave(filename = paste0(graph_save_location,"/",name_file,"_all_pheno.",graph_type),
                      plot = all_pheno_graph, device = graph_type, dpi = 300, width = 9,height = 6,units = "in")
    }
  }
}
#' Makes graphs from converted RDS files.
#'
#' @param x Group name
#' @param y Result file
#' @param save_folder Location of folder to save results
#' @param PheWAS_ID_filter List of PheWAS_IDs to filter tables and graph production.
#' @param plink_results Logical value if results are from Plink.
#' @param save_raw_plink Logical value indicating if wanting to save the raw plink results.
#' @param analysis_name name of analysis used in non-GRS data only.
#' @param SNP_list The list of SNPs used in extract_SNPs script (if used).
#' @param group_name_filter Name of group(s) to filter results by.
#' @param no_save_table_all_results logical if True then do not save all results.
#' @param SNP_filter List of SNPs to filter the graphs and tables production to.
#' @param MAC Minor allele count filter.
#' @param PheWAS_manifest PheWAS manifest data frame.
#' @param updated_manifest Altered PheWAS manifest to create all possible phenotypes.
#' @param sex_split Logical, if True three graphs are produced per input two sex stratified graphs and one graph with the original phenotype.
#' @param save_table_per_group_name Logical if True saves tables per group.
#' @param save_table_per_snp Logical if True saves tables per SNP.
#' @param R_association_graph Logical if True indicates use different graph input and data is from R-association rather than Plink.
#' @param graph_choice Value that indicates whether to produce one or two graphs and which ones.
#' @param R_association_results Logical if True indicates data is from R-association rather than Plink.
#' @param per_group_name_graph Logical if true produce graphs grouping by group.
#' @param per_snp_graph Logical if true produce graphs grouping by SNP
#' @param max_pheno_tests the maximum number of tests that are being run across grouping variable. Must be at => the number of associations.
#' @param FDR_figure value of significance for FDR.
#' @param max_FDR_graph value for lowest fdr value, converts 0 into this number.
#' @param MAC_figure Minor allele count filter numeric.
#' @param MAC_cases_N Minor allele count filter in cases
#' @param MAC_control_N Minor allele count filter in controls
#' @param PheWAS_label_filter list of PheWAS_IDs to label within the all_pheno graph.
#' @param max_overlap maximum number of overlaps for labelled phenotypes in the all_pheno graphs
#' @param graph_type save format of the graphs any input readable from ggsave is accepted
#' @param label_size size of the text for labelled phenotypes.
#' @param order_groups_alphabetically T or F to order groups on graph alphabetically rather than default lowest FDR.
#' @param order_phenotypes_alphabetically T or F to order phenotypes on graphs alphabetically rather than default lowest FDR.
#' @param save_all_graphs T or F to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return Saved tables and graphs, graphs created through calling additional functions.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
Deep_PheWAS_graphs_tables <- function(x,y,
                                      save_folder,
                                      PheWAS_ID_filter,
                                      plink_results,
                                      save_raw_plink,
                                      analysis_name,
                                      SNP_list,
                                      group_name_filter,
                                      no_save_table_all_results,
                                      SNP_filter,
                                      MAC,
                                      PheWAS_manifest,
                                      updated_manifest,
                                      sex_split,
                                      save_table_per_group_name,
                                      save_table_per_snp,
                                      R_association_graph,
                                      graph_choice,
                                      R_association_results,
                                      per_group_name_graph,
                                      per_snp_graph,
                                      max_pheno_tests,
                                      FDR_figure,
                                      max_FDR_graph,
                                      MAC_figure,
                                      MAC_cases_N,
                                      MAC_control_N,
                                      PheWAS_label_filter,
                                      max_overlap,
                                      graph_type,
                                      label_size,
                                      order_groups_alphabetically,
                                      order_phenotypes_alphabetically,
                                      save_all_graphs){

  . <- NULL
  # define save folder
  save_root <- save_folder

  dir.create(save_root,recursive = T)

  results <- y

  if(is.null(PheWAS_ID_filter)){
    PheWAS_ID_filter <- unique(results$PheWAS_ID)
  } else {
    PheWAS_IDs <- data.table::fread(PheWAS_ID_filter, header = F) %>%
      dplyr::pull(1)
    PheWAS_ID_filter <- unique(PheWAS_IDs)
  }

  results_PheWAS_ID_filter <- results %>%
    dplyr::filter(.data$PheWAS_ID %in% PheWAS_ID_filter)

  if(plink_results){
    if(save_raw_plink){
      save_name_plink <- paste0(analysis_name,"_",x,"_plink_results_raw.csv")
      data.table::fwrite(results,paste0(save_root,"/",save_name_plink))
    }

    if(all(!is.na(match(c("A1_CASE_CT","T_STAT"),colnames(results_PheWAS_ID_filter))))){

      main_table <- results_PheWAS_ID_filter %>%
        dplyr::left_join(SNP_list) %>%
        dplyr::filter(.data$ID %in% SNP_list$ID,
                      .data$rsid %in% SNP_filter,
                      .data$group_name %in% group_name_filter) %>%
        dplyr::mutate(dplyr::across(c(.data$POS,.data$A1_CT,.data$ALLELE_CT,.data$A1_CASE_CT,.data$A1_CTRL_CT,
                                      .data$A1_FREQ,.data$A1_CASE_FREQ,.data$A1_CTRL_FREQ,.data$MACH_R2,.data$OBS_CT,
                                      .data$OR,.data$`LOG(OR)_SE`,.data$L95,.data$U95,.data$Z_STAT,.data$P,.data$BETA,
                                      .data$SE,.data$T_STAT), as.numeric),
                      Beta=.data$BETA,
                      group=paste0(x),
                      minor_allele=ifelse(.data$A1_FREQ<0.5,.data$A1,.data$REF),
                      MAF=ifelse(.data$A1_FREQ<0.5,.data$A1_FREQ,1-.data$A1_FREQ),
                      SE=ifelse(is.na(.data$SE),.data$`LOG(OR)_SE`,.data$SE),
                      alter_direction=dplyr::case_when(.data$coded_allele==.data$A1 ~ 0,
                                                       .data$coded_allele!=.data$A1 ~ 1),
                      Beta=ifelse(is.na(.data$Beta),NA,ifelse(.data$alter_direction==1,-1*.data$Beta,.data$Beta)),
                      OR=ifelse(is.na(.data$OR),NA,ifelse(.data$alter_direction==1,1/.data$OR,.data$OR)),
                      N_L95=ifelse(is.na(.data$OR),ifelse(.data$alter_direction==1,-1*.data$U95,.data$L95),
                                   ifelse(.data$alter_direction==1,1/.data$U95,.data$L95)),
                      N_U95=ifelse(is.na(.data$OR),ifelse(.data$alter_direction==1,-1*.data$L95,.data$U95),
                                   ifelse(.data$alter_direction==1,1/.data$L95,.data$U95)),
                      effect_direction=ifelse(is.na(.data$OR),ifelse(is.na(.data$Beta),NA,ifelse(.data$Beta>0,"positive","negative")),
                                              ifelse(.data$OR>1,"positive","negative")),
                      Z_STAT=ifelse(.data$alter_direction==1,-1*.data$Z_STAT,.data$Z_STAT),
                      T_STAT=ifelse(.data$alter_direction==1,-1*.data$T_STAT,.data$T_STAT),
                      Z_T_STAT=ifelse(is.na(.data$Z_STAT),.data$T_STAT,.data$Z_STAT),
                      Error_flag=ifelse(.data$ERRCODE!=".",.data$ERRCODE,NA),
                      ID=as.factor(.data$ID),
                      MAC=ifelse(.data$A1_FREQ<0.5,.data$A1_CT,.data$ALLELE_CT-.data$A1_CT),
                      MAC_cases=ifelse(.data$A1_FREQ<0.5,.data$A1_CASE_CT,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*(1-.data$A1_CASE_FREQ)),
                      MAC_controls=ifelse(.data$A1_FREQ<0.5,.data$A1_CTRL_CT,(.data$A1_CTRL_CT/.data$A1_CTRL_FREQ)*(1-.data$A1_CTRL_FREQ)),
                      expected_MAC_cases=ifelse(.data$A1_FREQ<0.5,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*.data$A1_CTRL_FREQ,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*(1-.data$A1_CTRL_FREQ)),
                      MAC_diff=.data$MAC_cases-.data$expected_MAC_cases,
                      keep=ifelse(!is.na(.data$Beta),2,ifelse(.data$expected_MAC_cases<=2&.data$MAC_cases>=3,1,ifelse(.data$expected_MAC_cases>=7&(.data$MAC_cases<=5&.data$MAC_cases>=1),1,0))),
                      ratio=1/((.data$A1_CASE_CT/.data$A1_CASE_FREQ)/(.data$A1_CTRL_CT/.data$A1_CTRL_FREQ)),
                      sex_pheno=sub(".*_", "", .data$PheWAS_ID),
                      join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female")) %>%
        dplyr::filter(.data$MAC >=MAC_figure,
                      dplyr::case_when(.data$keep==0 ~ .data$MAC_cases >= MAC_cases_N & .data$MAC_controls >= MAC_control_N,
                                       .data$keep==1 ~ .data$MAC_cases>=1 & .data$MAC_controls >= MAC_control_N,
                                       .data$keep==2 ~ .data$MAC >= MAC_figure)) %>%
        dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
        dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc)) %>%
        dplyr::select(.data$group,collective_name=.data$group_name,.data$PheWAS_ID,.data$category,description=.data$phenotype,N_ID=.data$OBS_CT,.data$rsid,.data$P,.data$OR,.data$Beta,L95=.data$N_L95,U95=.data$N_U95,.data$coded_allele,.data$non_coded_allele,.data$minor_allele,.data$MAF,.data$MAC,.data$MAC_cases,.data$MAC_controls,chromosome=.data$`#CHROM`,position=.data$POS,.data$Z_T_STAT,.data$SE,.data$effect_direction,.data$category,phenotype_group=.data$pheno_group,phenoytpe_group_narrow=.data$group_narrow,.data$short_desc,Info_score=.data$MACH_R2,firth=.data$`FIRTH?`,.data$TEST,.data$Error_flag,.data$ID,.data$graph_save_name)

    } else if(all(!is.na(match("T_STAT",colnames(results_PheWAS_ID_filter))) &&  is.na(match("A1_CASE_CT",colnames(results_PheWAS_ID_filter))))){

      main_table <- results_PheWAS_ID_filter %>%
        dplyr::left_join(SNP_list) %>%
        dplyr::filter(.data$ID %in% SNP_list$ID,
                      .data$rsid %in% SNP_filter,
                      .data$group_name %in% group_name_filter) %>%
        dplyr::mutate(dplyr::across(c(.data$POS,.data$A1_CT,.data$ALLELE_CT,.data$A1_FREQ,.data$MACH_R2,.data$OBS_CT,
                                      .data$L95,.data$U95,.data$P,.data$BETA,
                                      .data$SE,.data$T_STAT), as.numeric),
                      Beta=.data$BETA,
                      group=paste0(x),
                      minor_allele=ifelse(.data$A1_FREQ<0.5,.data$A1,.data$REF),
                      MAF=ifelse(.data$A1_FREQ<0.5,.data$A1_FREQ,1-.data$A1_FREQ),
                      alter_direction=dplyr::case_when(.data$coded_allele==.data$A1 ~ 0,
                                                       .data$coded_allele!=.data$A1 ~ 1),
                      Beta=ifelse(is.na(.data$Beta),NA,ifelse(.data$alter_direction==1,-1*.data$Beta,.data$Beta)),
                      N_L95=ifelse(.data$alter_direction==1,-1*.data$U95,.data$L95),
                      N_U95=ifelse(.data$alter_direction==1,-1*.data$L95,.data$U95),
                      effect_direction=ifelse(is.na(.data$Beta),NA,ifelse(.data$Beta>0,"positive","negative")),
                      Z_T_STAT=ifelse(.data$alter_direction==1,-1*.data$T_STAT,.data$T_STAT),
                      Error_flag=ifelse(.data$ERRCODE!=".",.data$ERRCODE,NA),
                      ID=as.factor(.data$ID),
                      MAC=ifelse(.data$A1_FREQ<0.5,.data$A1_CT,.data$ALLELE_CT-.data$A1_CT),
                      keep=ifelse(!is.na(.data$Beta),2,0),
                      sex_pheno=sub(".*_", "", .data$PheWAS_ID),
                      join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female")) %>%
        dplyr::filter(.data$MAC >=MAC_figure,
                      dplyr::case_when(.data$keep==2 ~ .data$MAC >= MAC_figure)) %>%
        dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
        dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc)) %>%
        dplyr::select(.data$group,collective_name=.data$group_name,.data$PheWAS_ID,.data$category,description=.data$phenotype,N_ID=.data$OBS_CT,.data$rsid,.data$P,.data$Beta,L95=.data$N_L95,U95=.data$N_U95,.data$coded_allele,.data$non_coded_allele,.data$minor_allele,.data$MAF,.data$MAC,chromosome=.data$`#CHROM`,position=.data$POS,.data$Z_T_STAT,.data$SE,.data$effect_direction,.data$category,phenotype_group=.data$pheno_group,phenoytpe_group_narrow=.data$group_narrow,.data$short_desc,Info_score=.data$MACH_R2,.data$TEST,.data$Error_flag,.data$ID,.data$graph_save_name)


    } else if(all(!is.na(match("A1_CASE_CT",colnames(results_PheWAS_ID_filter))) &&  is.na(match("T_STAT",colnames(results_PheWAS_ID_filter))))){
      main_table <- results_PheWAS_ID_filter %>%
        dplyr::left_join(SNP_list) %>%
        dplyr::filter(.data$ID %in% SNP_list$ID,
                      .data$rsid %in% SNP_filter,
                      .data$group_name %in% group_name_filter) %>%
        dplyr::mutate(dplyr::across(c(.data$POS,.data$A1_CT,.data$ALLELE_CT,.data$A1_CASE_CT,.data$A1_CTRL_CT,
                                      .data$A1_FREQ,.data$A1_CASE_FREQ,.data$A1_CTRL_FREQ,.data$MACH_R2,.data$OBS_CT,
                                      .data$OR,.data$`LOG(OR)_SE`,.data$L95,.data$U95,.data$Z_STAT,.data$P), as.numeric),
                      group=paste0(x),
                      minor_allele=ifelse(.data$A1_FREQ<0.5,.data$A1,.data$REF),
                      MAF=ifelse(.data$A1_FREQ<0.5,.data$A1_FREQ,1-.data$A1_FREQ),
                      SE=.data$`LOG(OR)_SE`,
                      alter_direction=dplyr::case_when(.data$coded_allele==.data$A1 ~ 0,
                                                       .data$coded_allele!=.data$A1 ~ 1),
                      OR=ifelse(.data$alter_direction==1,1/.data$OR,.data$OR),
                      N_L95=ifelse(.data$alter_direction==1,1/.data$U95,.data$L95),
                      N_U95=ifelse(.data$alter_direction==1,1/.data$L95,.data$U95),
                      effect_direction=ifelse(.data$OR>1,"positive","negative"),
                      Z_T_STAT=ifelse(.data$alter_direction==1,-1*.data$Z_STAT,.data$Z_STAT),
                      Error_flag=ifelse(.data$ERRCODE!=".",.data$ERRCODE,NA),
                      ID=as.factor(.data$ID),
                      MAC=ifelse(.data$A1_FREQ<0.5,.data$A1_CT,.data$ALLELE_CT-.data$A1_CT),
                      MAC_cases=ifelse(.data$A1_FREQ<0.5,.data$A1_CASE_CT,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*(1-.data$A1_CASE_FREQ)),
                      MAC_controls=ifelse(.data$A1_FREQ<0.5,.data$A1_CTRL_CT,(.data$A1_CTRL_CT/.data$A1_CTRL_FREQ)*(1-.data$A1_CTRL_FREQ)),
                      expected_MAC_cases=ifelse(.data$A1_FREQ<0.5,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*.data$A1_CTRL_FREQ,(.data$A1_CASE_CT/.data$A1_CASE_FREQ)*(1-.data$A1_CTRL_FREQ)),
                      MAC_diff=.data$MAC_cases-.data$expected_MAC_cases,
                      keep=ifelse(.data$expected_MAC_cases<=2&.data$MAC_cases>=3,1,ifelse(.data$expected_MAC_cases>=7&(.data$MAC_cases<=5&.data$MAC_cases>=1),1,0)),
                      ratio=1/((.data$A1_CASE_CT/.data$A1_CASE_FREQ)/(.data$A1_CTRL_CT/.data$A1_CTRL_FREQ)),
                      sex_pheno=sub(".*_", "", .data$PheWAS_ID),
                      join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female")) %>%
        dplyr::filter(.data$MAC >=MAC_figure,
                      dplyr::case_when(.data$keep==0 ~ .data$MAC_cases >= MAC_cases_N & .data$MAC_controls >= MAC_control_N,
                                       .data$keep==1 ~ .data$MAC_cases>=1 & .data$MAC_controls >= MAC_control_N,
                                       .data$keep==2 ~ .data$MAC >= MAC_figure)) %>%
        dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
        dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc)) %>%
        dplyr::select(.data$group,collective_name=.data$group_name,.data$PheWAS_ID,.data$category,description=.data$phenotype,N_ID=.data$OBS_CT,.data$rsid,.data$P,.data$OR,L95=.data$N_L95,U95=.data$N_U95,.data$coded_allele,.data$non_coded_allele,.data$minor_allele,.data$MAF,.data$MAC,.data$MAC_cases,.data$MAC_controls,chromosome=.data$`#CHROM`,position=.data$POS,.data$Z_T_STAT,.data$SE,.data$effect_direction,.data$category,phenotype_group=.data$pheno_group,phenoytpe_group_narrow=.data$group_narrow,.data$short_desc,Info_score=.data$MACH_R2,firth=.data$`FIRTH?`,.data$TEST,.data$Error_flag,.data$ID,.data$graph_save_name)

    } else {
      rlang::abort(paste0("Results loaded do not have the correct columns and are missing at minimum 'T_STAT' or 'A1_CASE_CT'"))
    }
    main_table_fdr_split <- main_table %>%
      dplyr::filter(.data$PheWAS_ID %in% updated_manifest$PheWAS_ID) %>%
      dplyr::group_split(.data$ID)
  }
  if(R_association_results){
    main_table <- results_PheWAS_ID_filter %>%
      dplyr::filter(.data$PheWAS_ID %in% updated_manifest$PheWAS_ID) %>%
      dplyr::mutate(dplyr::across(c(.data$P,.data$OR,.data$Beta,.data$SE,.data$L95,.data$U95), as.numeric)) %>%
      dplyr::select(-.data$phenotype,-.data$phenotype_group,-.data$group_narrow,-.data$short_desc) %>%
      dplyr::left_join(updated_manifest,by="PheWAS_ID")%>%
      dplyr::filter(.data$PheWAS_ID %in% updated_manifest$PheWAS_ID) %>%
      dplyr::select(.data$name,.data$PheWAS_ID,.data$category,.data$phenotype,.data$P,.data$OR,.data$Beta,.data$L95,.data$U95,.data$SE,.data$phenotype_group,.data$group_narrow,.data$short_desc,.data$effect_direction,.data$group,.data$name_group)

    main_table_fdr_split <- tibble::as_tibble(main_table) %>%
      dplyr::group_split(.data$name)
  }
  if(length(main_table_fdr_split)==0){
    return()
  }
  main_table_fdr <- mapply(fdr_calc,main_table_fdr_split,MoreArgs = list(max_pheno_tests=max_pheno_tests),SIMPLIFY = F) %>%
    data.table::rbindlist(.) %>%
    dplyr::relocate(.data$FDR,.before = .data$P) %>%
    dplyr::arrange(.data$FDR) %>%
    dplyr::mutate(sex_pheno_identifier="")
  if(sex_split){
    main_table_fdr <- main_table_fdr %>%
      dplyr::mutate(male=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male$"),1,0),
                    female=ifelse(stringr::str_detect(.data$PheWAS_ID,"_female$"),1,0),
                    combined=ifelse(.data$male==1 | .data$female ==1,0,1),
                    sex_pheno_identifier=dplyr::case_when(.data$male==1 ~ "male",
                                                          .data$female==1 ~ "female",
                                                          .data$combined==1 ~ "")) %>%
      dplyr::select(-.data$female,-.data$male,-.data$combined)
    message(paste("loc_3"))
  }

  main_table_fdr_split <- main_table_fdr %>%
    dplyr::group_split(.data$sex_pheno_identifier)

  if(isFALSE(no_save_table_all_results)) {
    data.table::fwrite(main_table_fdr,paste0(save_root,"/",analysis_name,"_",x,"_filtered_results.csv"))
  }

  if(plink_results){
    if(save_table_per_group_name) {
      per_group_name <- main_table_fdr %>%
        dplyr::group_by(.data$collective_name,.data$PheWAS_ID) %>%
        dplyr::slice_min(.data$FDR)
      per_group_name_split <- per_group_name %>%
        dplyr::ungroup() %>%
        dplyr::group_split(.data$collective_name)

      grouping <- x
      per_group_name_folder <- paste0(save_root,"/",analysis_name,"_",x,"_results_per_collective_name/")
      dir.create(per_group_name_folder)
      data.table::fwrite(per_group_name,paste0(save_root,x,"_results_per_collective_name_all_",analysis_name,".csv"))
      lapply(per_group_name_split,function(x) data.table::fwrite(x,paste0(per_group_name_folder,unique(x$collective_name),"_",grouping,"_results.csv")))
    }
    if(save_table_per_snp){
      per_snp_tables <- main_table_fdr %>%
        dplyr::group_split(.data$ID)
      grouping <- x
      per_SNP_folder <- paste0(save_root,"/",analysis_name,"_",x,"_results_per_SNP/")
      dir.create(per_SNP_folder)
      lapply(per_snp_tables,function(x) data.table::fwrite(x,paste0(per_SNP_folder,unique(x$graph_save_name),"_",grouping,"_results.csv")))
    }
    if(per_group_name_graph) {

      mapply(per_group_function,main_table_fdr_split,MoreArgs = list(y=x,save_root=save_root,analysis_name=analysis_name,graph_choice=graph_choice,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,PheWAS_label_filter=PheWAS_label_filter,max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                                                     order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                                     save_all_graphs=save_all_graphs))
    }
    if(per_snp_graph) {

      mapply(per_snp_function,main_table_fdr_split,MoreArgs = list(y=x,save_root=save_root,analysis_name=analysis_name,graph_choice=graph_choice,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,PheWAS_label_filter=PheWAS_label_filter,max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                                                   order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                                   save_all_graphs=save_all_graphs))
    }
  }

  if(R_association_graph) {
    mapply(R_association_function,main_table_fdr_split,MoreArgs = list(save_root=save_root,analysis_name=analysis_name,
                                                                       graph_choice=graph_choice,FDR_figure=FDR_figure,max_FDR_graph=max_FDR_graph,PheWAS_label_filter=PheWAS_label_filter,
                                                                       max_overlap=max_overlap,graph_type=graph_type,label_size=label_size,order_groups_alphabetically=order_groups_alphabetically,
                                                                       order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                                       save_all_graphs=save_all_graphs))
  }

}


#' Makes graphs and tables from RDS results files.
#'
#'@param results_file Full file path of the results file RDS R list object.
#'@param analysis_name Name for the analysis.
#'@param plink_results Select if results are from 03a_PLINK_association_testing.R
#'@param SNP_list Full path of the SNP_list file used in 02_extracting_snps.R.
#'@param R_association_results Select if results are from 03b_R_association_testing.R.
#'@param save_folder File path of the folder to which the output will be saved.
#'@param group_filter Comma-separated text input, used to filter the group to which the table and graph functions are applied. Inputted groups are the ones that are retained for analysis, group here refers to the grouping variable used to subset the analysis classically ancestry.
#'@param PheWAS_ID_filter Full file path to a plain txt file containing single column NO header containing full PheWAS_ID of phenotypes that will be included Use if wanting to apply the table and/or graphing functions a subset of phenotypes.
#'@param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#'@param max_pheno Manual override for inputting maximum phenotypes analysed. Used for calculating FDR. The default used the largest number of associations in the results_file of any grouping.
#'@param sig_FDR Value of FDR for which associations are reported as significant. Will alter the line of significance in both graph types and the reported phenotypes in the significant_pheno graphs. Default: 0.01.
#'@param no_save_table_all_results Specify whether a table containing all the results should be saved (FALSE) or not (TRUE). Default is FALSE.
#'@param no_graph_all Specify whether a figure of all the results should be produced (FALSE) or not (TRUE). Default is FALSE.
#'@param no_graph_sig Specify whether a figure of the significant results should be produced (FALSE) or not (TRUE). Default is FALSE.
#'@param max_FDR_graph Maximum value for the y-axis to be used when the FDR for a particular association is exceptionally small. Default is 300.
#'@param sex_split Specify whether the figure should be stratified by sex (TRUE) or not (FALSE). If TRUE, this will produce separate figures for males, females and combined. Default is FALSE. Does not create split tables.
#'@param SNP_filter Full path of a file containing SNP IDs, used to filter the output to a given subset of SNPs. Use if wanting to apply the table and/or graphing functions a subset of results. Only works for plink_results.
#'@param group_name_filter Full path of the file containing a list of group_names that correspond with those defined in SNP_list. This can be used if the table and/or graphing functions should be applied to a subset of SNPs, such as only those in a given credible set. This option can only be used with association results generated by 03a_PLINK_association_testing.R.
#'@param save_raw_plink  Specify whether the unfiltered association results from the PLINK analysis should be saved (TRUE) or not (FALSE). Default is FALSE.
#'@param MAC MAC (minor allele count) filter applied to all associations, only applicable in results from03a_PLINK_association_testing.R Default is 20.
#'@param MAC_case Minor allele count among cases for filtering the association results. only applicable in results from 03a_PLINK_association_testing.R and only in binary phenotypes. Default is 5.
#'@param MAC_control Minor allele count among controls for filtering the association results. Only applicable in results from 03a_PLINK_association_testing.R and only in binary phenotypes. Default: 10.
#'@param per_group_name_graph Select if wanting to produce graphs per-group_name. This is used when looking to report the most significant finding across several SNPs for a single construct, potentially and gene or a sentinel SNP with a credible set. It is a column in the SNP_list file. Default is FALSE.
#'@param per_snp_graph Specify whether a figure is produced for every SNP provided (TRUE) or not (FALSE). Default is FALSE.
#'@param save_table_per_group_name Select if wanting to produce tables per-group_name. This is used when looking to report the most significant finding across several SNPs for a single construct, potentially and gene or a sentinel SNP with a credible set. It is a column in the SNP_list file. Default is FALSE.
#'@param save_table_per_snp Specify whether a results table should be generated for every SNP provided (TRUE) or not (FALSE). Default is FALSE. Will be saved in a created folder named /analysis_name_group_per_SNP_tables. Example if the group was groupA and analysis_name top_SNPs the folder would be /top_SNPs_groupA_per_SNP_tables.
#'@param R_association_graph Specify whether a figure of the results from the association analysis from 03b_R_association_testing.R should be produced (TRUE) or not (FALSE). Default is FALSE.
#'@param PheWAS_ID_label_filter Full file path to plain text file containing single headed column of PheWAS_IDs. Only these PheWAS_IDs will be labelled within any graphical output designed to be used primarily when trying to edit a single graph as the filter will apply to all graphs being created.
#'@param max_overlap_labels Number, represents maximum overlaps for labelling of phenotypes in the all_pheno graph, lowering the number has the effect of reducing the total number of phenotypes labelled. Default: 20
#'@param graph_file_save Allows user to specify the file format of the graphs, for example pdf or png. Default: png
#'@param label_text_size Number, represents the text size of the labelled phenotypes in all_pheno graph. Default: 2
#'@param order_groups_alphabetically Specify whether the groups in the graphs should be ordered by lowest FDR (FALSE) or alphabetically (TRUE).
#'@param order_phenotypes_alphabetically Specify whether the phenotypes within each group are ordered by the lowest FDR (FALSE) or alphabetically (TRUE).
#'@param save_all_graphs Specify whether to always save every graph with or without a significant result (TRUE) or to only save when at least one association is significant (FALSE).
#' @return Tables and graphs of the association results.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
graphs_tables_DeepPheWAS <- function(results_file,
                                     analysis_name,
                                     plink_results,
                                     SNP_list,
                                     R_association_results,
                                     save_folder,
                                     group_filter,
                                     PheWAS_ID_filter,
                                     PheWAS_manifest_overide,
                                     max_pheno,
                                     sig_FDR,
                                     no_save_table_all_results,
                                     no_graph_all,
                                     no_graph_sig,
                                     max_FDR_graph,
                                     sex_split,
                                     SNP_filter,
                                     group_name_filter,
                                     save_raw_plink,
                                     MAC,
                                     MAC_case,
                                     MAC_control,
                                     per_group_name_graph,
                                     per_snp_graph,
                                     save_table_per_group_name,
                                     save_table_per_snp,
                                     R_association_graph,
                                     PheWAS_ID_label_filter,
                                     max_overlap_labels,
                                     graph_file_save,
                                     label_text_size,
                                     order_groups_alphabetically,
                                     order_phenotypes_alphabetically,
                                     save_all_graphs){
  if(is.null(PheWAS_manifest_overide)){
    PheWAS_manifest <- data.table::fread(system.file("extdata","PheWAS_manifest.csv.gz", package = "DeepPheWAS"))
  } else {
    if(!file.exists(PheWAS_manifest_overide)){
      rlang::abort(paste0("'PheWAS_manifest_overide' must be a file"))
    }
    PheWAS_manifest <- data.table::fread(PheWAS_manifest_overide)
    PheWAS_manifest_overide_colname <- c("PheWAS_ID","broad_catagory","category","phenotype","range_ID","sex","field_code","QC_flag_ID","QC_filter_values","first_last_max_min_value","date_code","age_code","case_code","control_code","exclude_case_control_crossover","included_in_analysis","quant_combination","primary_care_code_list","limits","lower_limit","upper_limit","transformation","analysis","age_col","pheno_group","short_desc","group_narrow","concept_name","concept_description","Notes")

    if(!all(tibble::has_name(PheWAS_manifest,PheWAS_manifest_overide_colname))){
      warning(paste0("'PheWAS_manifest_overide' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(PheWAS_manifest_overide_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(PheWAS_manifest), collapse=","),
              paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(PheWAS_manifest),PheWAS_manifest_overide_colname), collapse=","))
    }
  }

  # create an updated manifests that includes sex-specific phenotypes, equates ot maximum possible phenotypes
  PheWAS_manifest_edit <- PheWAS_manifest %>%
    dplyr::select(.data$PheWAS_ID,.data$sex,.data$phenotype,.data$category,.data$pheno_group,.data$group_narrow,.data$short_desc,.data$included_in_analysis) %>%
    dplyr::filter(.data$included_in_analysis==1)
  sex_names <- c("Male","Female")
  PheWAS_manifest_male <- PheWAS_manifest_edit %>%
    dplyr::filter(!.data$sex %in% sex_names) %>%
    dplyr::mutate(PheWAS_ID=paste0(.data$PheWAS_ID,"_male"),
                  short_desc=paste(.data$short_desc," (male)"),
                  phenotype=paste0("Male stratified variant of ",.data$phenotype))
  PheWAS_manifest_female <- PheWAS_manifest_edit %>%
    dplyr::filter(!.data$sex %in% sex_names) %>%
    dplyr::mutate(PheWAS_ID=paste0(.data$PheWAS_ID,"_female"),
                  short_desc=paste(.data$short_desc," (female)"),
                  phenotype=paste0("Female stratified variant of ",.data$phenotype))
  PheWAS_manifest_AOO <-  PheWAS_manifest_edit %>%
    dplyr::mutate(PheWAS_ID=paste0(.data$PheWAS_ID,"_age_of_onset"),
                  short_desc=paste(.data$short_desc," (age of onset)"),
                  phenotype=paste0("Age of onset of ",.data$phenotype))

  updated_manifest <- PheWAS_manifest_edit %>%
    dplyr::bind_rows(PheWAS_manifest_male,PheWAS_manifest_female,PheWAS_manifest_AOO) %>%
    dplyr::select(.data$PheWAS_ID,.data$phenotype,.data$category,phenotype_group=.data$pheno_group,.data$group_narrow,.data$short_desc)
  #read in all results
  all_results <- purrr::compact(readRDS(results_file))
  # group_filter
  if(!is.null(group_filter)) {
    group_filter <- unlist(strsplit(group_filter,","))
  } else {
    group_filter <- names(all_results)
  }
  all_results <- all_results[group_filter]
  group_names <- names(all_results)
  # edit to plink results
  if(plink_results){
    # SNP_list edit only for plink list
    SNP_list <- data.table::fread(SNP_list) %>%
      dplyr::mutate(group=.data$group_name,
                    coded_number=match(stringr::str_sub(.data$coded_allele,1,1),LETTERS),
                    non_coded_number=match(stringr::str_sub(.data$non_coded_allele,1,1),LETTERS),
                    multi=dplyr::case_when(nchar(.data$coded_allele) > 1 & nchar(.data$non_coded_allele)== 1 ~ 1,
                                           nchar(.data$non_coded_allele) > 1 & nchar(.data$coded_allele)== 1 ~ 2,
                                           nchar(.data$coded_allele) > 1 & nchar(.data$non_coded_allele) <1 ~ 3),
                    allele_order=ifelse(.data$non_coded_number>.data$coded_number,paste0(.data$coded_allele,"_",.data$non_coded_allele),
                                        ifelse(.data$coded_number>.data$non_coded_number,paste0(.data$non_coded_allele,"_",.data$coded_allele),
                                               ifelse(.data$coded_number==.data$non_coded_number,
                                                      ifelse(.data$multi==1,paste0(.data$non_coded_allele,"_",.data$coded_allele),
                                                             ifelse(.data$multi==2,paste0(.data$coded_allele,"_",.data$non_coded_allele),
                                                                    ifelse(nchar(.data$coded_allele)>nchar(.data$non_coded_allele),paste0(.data$non_coded_allele,"_",.data$coded_allele),paste0(.data$coded_allele,"_",.data$non_coded_allele)))),NA))),
                    ID=paste0(.data$rsid,"_",.data$allele_order)) %>%
      dplyr::select(.data$chromosome,.data$rsid,.data$group_name,.data$coded_allele,.data$non_coded_allele,.data$graph_save_name,.data$ID)
    # create SNP_filter
    if(!is.null(SNP_filter)) {
      SNP_filter <- data.table::fread(SNP_filter) %>%
        dplyr::pull(.data$rsid)
    } else {
      SNP_filter <- SNP_list %>%
        dplyr::pull(.data$rsid)
    }
    # group_name filter
    if(!is.null(group_name_filter)) {
      group_name_filter <- data.table::fread(group_name_filter) %>%
        dplyr::pull(.data$group_name)
    } else {
      group_name_filter <- SNP_list %>%
        dplyr::pull(.data$group_name)
    }
  }
  # inputs into graphs and tables limits to only those that should have been analysed, but does not limit to those selected for graphs. Or is manually selected
  if(is.null(max_pheno)) {
    max_pheno_tests <- max(unlist(lapply(all_results, function(x) nrow(data.frame(ids=unique(x$PheWAS_ID)) %>% dplyr::filter(.data$ids %in% updated_manifest$PheWAS_ID)))))
  } else {
    max_pheno_tests <- as.numeric(max_pheno)
  }
  if(is.na(as.numeric(sig_FDR))){
    rlang::abort(paste0("'sig_FDR' must be a numeral"))
  }
  FDR_figure <- as.numeric(sig_FDR)
  if(is.na(as.numeric(MAC))){
    rlang::abort(paste0("'MAC' must be a numeral"))
  }
  MAC_figure <- as.numeric(MAC)
  if(is.na(as.numeric(MAC_case))){
    rlang::abort(paste0("'MAC_case' must be a numeral"))
  }
  MAC_cases_N <- as.numeric(MAC_case)
  if(is.na(as.numeric(MAC_control))){
    rlang::abort(paste0("'MAC_control' must be a numeral"))
  }
  MAC_control_N <- as.numeric(MAC_control)
  if(is.na(as.numeric(max_FDR_graph))){
    rlang::abort(paste0("'max_FDR_graph' must be a numeral"))
  }
  max_FDR_graph_parse <- as.numeric(max_FDR_graph)
  max_FDR_graph <- 10^(-max_FDR_graph_parse)
  if(!is.null(PheWAS_ID_label_filter)){
    if(!file.exists(PheWAS_ID_label_filter)){
      rlang::abort(paste0("'PheWAS_ID_label_filter' must be a file"))
    }
    PheWAS_label_filter <- data.table::fread(PheWAS_ID_label_filter) %>%
      dplyr::pull(1)
  } else {
    PheWAS_label_filter <- NULL
  }
  if(is.na(as.numeric(max_overlap_labels))){
    rlang::abort(paste0("'max_overlap_labels' must be a numeral"))
  }
  max_overlap <- as.numeric(max_overlap_labels)
  graph_type <- graph_file_save
  if(is.na(as.numeric(label_text_size))){
    rlang::abort(paste0("'label_text_size' must be a numeral"))
  }
  label_size <- as.numeric(label_text_size)

  # graph choices
  if(no_graph_all & isFALSE(no_graph_sig)) {
    graph_choice <- "sig_only"
  } else if(isFALSE(no_graph_all) & no_graph_sig) {
    graph_choice <- "all_pheno"
  } else {
    graph_choice <- "both"
  }

  mapply(Deep_PheWAS_graphs_tables,group_names,all_results,MoreArgs = list(save_folder=save_folder,
                                                                           PheWAS_ID_filter=PheWAS_ID_filter,
                                                                           plink_results=plink_results,
                                                                           save_raw_plink=save_raw_plink,
                                                                           analysis_name=analysis_name,
                                                                           SNP_list=SNP_list,
                                                                           group_name_filter=group_name_filter,
                                                                           no_save_table_all_results=no_save_table_all_results,
                                                                           SNP_filter=SNP_filter,
                                                                           MAC=MAC,
                                                                           PheWAS_manifest=PheWAS_manifest,
                                                                           updated_manifest=updated_manifest,
                                                                           sex_split=sex_split,
                                                                           save_table_per_group_name=save_table_per_group_name,
                                                                           save_table_per_snp=save_table_per_snp,
                                                                           R_association_graph=R_association_graph,
                                                                           graph_choice=graph_choice,
                                                                           R_association_results=R_association_results,
                                                                           per_group_name_graph=per_group_name_graph,
                                                                           per_snp_graph=per_snp_graph,
                                                                           max_pheno_tests=max_pheno_tests,
                                                                           FDR_figure=FDR_figure,
                                                                           max_FDR_graph=max_FDR_graph,
                                                                           MAC_figure=MAC_figure,
                                                                           MAC_cases_N=MAC_cases_N,
                                                                           MAC_control_N=MAC_control_N,
                                                                           PheWAS_label_filter=PheWAS_label_filter,
                                                                           max_overlap=max_overlap,
                                                                           graph_type=graph_type,
                                                                           label_size=label_size,
                                                                           order_groups_alphabetically=order_groups_alphabetically,
                                                                           order_phenotypes_alphabetically=order_phenotypes_alphabetically,
                                                                           save_all_graphs=save_all_graphs))
}
