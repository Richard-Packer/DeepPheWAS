#' Combines PLINK results from the split_groups and split_analysis.
#'
#' @param x Phenotype results file.
#' @param y full file location for save file minus specific file name
#' @param remove_string name of the analysis
#' @param analysis_name Name for the analysis to be used in the naming of the output files.
#' @return list of results
#' @importFrom magrittr %>%
plink_file_edit_join <- function(x,y,remove_string,analysis_name) {
  name <- data.frame(stringr::str_remove_all(y, paste(remove_string, collapse = "|"))) %>%
    dplyr::rename(PheWAS_ID=1)
  per_result <- fread(x) %>%
    bind_cols(name) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.character),
                  table_save=analysis_name)
  return(per_result)
}
#' Combines PLINK results from the split_groups and split_analysis into a single file.
#'
#' @param split_plink_results_folder Full path of the folder containing all of the PLINK results from the split_group, split_analysis option in 03a_PLINK_association_testing.R.
#' @param results_RDS_file Full path of the results of the analysis from the 03a_PLINK_association_testing.R script. Is an R object.
#' @param group_name Name of the group that was split that the results represent.
#' @param analysis_name Name for the analysis to be used in the naming of the output files, use the same as in 03a_PLINK_association_analysis.
#' @return an updated results RDS file that includes the now combined group results.
#' @importFrom magrittr %>%
#' @export
combine_split <- function(split_plink_results_folder,
                          results_RDS_file,
                          group_name,
                          analysis_name){
  . <- NULL
  # combine results
  results_files <-  list.files(path = split_plink_results_folder, pattern = paste0(group_name,"(.*?)hybrid|",group_name,"(.*?)linear"))
  remove_string <- c(".glm.logistic.hybrid", paste0(group_name,"."), ".glm.linear")
  results_file_location <- list.files(path = split_plink_results_folder, pattern = paste0(group_name,"(.*?)hybrid|",group_name,"(.*?)linear"),full.names = T)

  # join the plink results into a single table
  plink_results <- mapply(plink_file_edit_join,
                          results_file_location,
                          results_files,
                          MoreArgs = list(remove_string=remove_string,
                                          analysis_name=analysis_name),
                          SIMPLIFY = F) %>%
    data.table::rbindlist(.,fill = T)
  group_name <- group_name
  plink_results_list <- list(plink_results)
  names(plink_results_list) <- group_name
  # now combine with original list of results
  existing_results <- readRDS(results_RDS_file)
  # remove null
  existing_results_trim <- purrr::compact(existing_results)
  # append results
  new_results <- append(existing_results_trim,plink_results_list)
  saveRDS(new_results,results_RDS_file)
}
