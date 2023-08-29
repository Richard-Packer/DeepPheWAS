#' Summarises the number of cases and controls per_phenotype inputed.
#'
#' @param x Phenotype
#' @param y group
#' @param z PheWAS_ID.
#' @param PheWAS_manifest manifest file
#' @return A dataframe of summary stats for phenotypes.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
sum_stats <- function(x,y,z,PheWAS_manifest) {

  name <- stringr::str_remove(z,"_male|_female")

  if(stringr::str_detect(name,"age_of_onset")){
    type <- "quant"
  } else {
    type <- PheWAS_manifest %>%
      dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
      dplyr::select(.data$analysis) %>%
      dplyr::pull()
  }

  if (type=="quant") {
    cases <- x %>%
      tidyr::drop_na(2)

    n_cases <- nrow(cases)
    n_control <- NA

    phecode_name <- z
    stats_name <- paste0("N_cases_",y)
    control_name <- paste0("N_control_",y)

    summary_sample_stats <- data.frame(phecode_name,n_cases,n_control) %>%
      dplyr::rename(PheWAS_ID=1,!!stats_name:=2,!!control_name:=3)
  } else {
    ## epi cases and controls
    cases <- x
    n_cases <- length(which(cases[,2]>0))
    n_control <- length(which(cases[,2]==0))

    phecode_name <- z
    stats_name <- paste0("N_cases_",y)
    control_name <- paste0("N_control_",y)

    summary_sample_stats <- data.frame(phecode_name,n_cases,n_control) %>%
      dplyr::rename(PheWAS_ID=1,!!stats_name:=2,!!control_name:=3)
  }
  return(summary_sample_stats)
}

#' Creates grouping ID's
#'
#' @param x Groups
#' @param groupings_df grouping data frame
#' @return Grouping ID dataframes as lists
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
group_id_creater <- function(x,groupings_df) {
  ID <- groupings_df %>%
    dplyr::filter(.data$group==x) %>%
    dplyr::pull(.data$eid)
  return(ID)
}
#' Removes related individuals from phenotypes.
#'
#' @param x Phenotypes extracted
#' @param y PheWAS_ID
#' @param call_rate_kinships call_rate_kinship file
#' @param PheWAS_manifest manifest file
#' @return A list phenotypes filtered for relatedness.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
relate_remove_fn <- function(x,y,call_rate_kinships,PheWAS_manifest) {

  name <- stringr::str_remove(y,"_male|_female")

  if(stringr::str_detect(name,"age_of_onset")){
    type <- "quant"
  } else {
    type <- PheWAS_manifest %>%
      dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
      dplyr::select(.data$analysis) %>%
      dplyr::pull()
  }

  pheno <- x %>%
    dplyr::rename(any_code=2)

  if (type=="quant") {
    rel<-call_rate_kinships

    pheno<-x %>%
      tidyr::drop_na()

    # Remove all rows from the relatedness file where (1) one or both IDs are not in the phenotype
    #(1)
    rel <- call_rate_kinships %>%
      dplyr::filter(.data$ID1 %in% pheno$eid & .data$ID2 %in% pheno$eid)

    # Create identifier for whether pair contains a duplicated ID or not and partition into two data frames (independent pairs and pairs containing duplicated IDs)
    AllID<-c(rel$ID1,rel$ID2)
    dupID<-AllID[which(duplicated(AllID))]
    rel$dup<-"NO"
    rel$dup[which(rel$ID1%in%dupID|rel$ID2%in%dupID)]<-"YES"

    indpairs<-rel[which(rel$dup=="NO"),]
    duppairs<-rel[which(rel$dup=="YES"),]

    # How many times does each duplicated ID arise?

    repeat_counter <- function(x) {

      duplicated_count_ID1 <- rel %>%
        dplyr::filter(.data$dup=="YES") %>%
        dplyr::select(ID=.data$ID1)

      duplicated_count_ID2 <- rel %>%
        dplyr::filter(.data$dup=="YES") %>%
        dplyr::select(ID=.data$ID2)

      duplicated_count_all <- duplicated_count_ID1 %>%
        dplyr::bind_rows(duplicated_count_ID2) %>%
        dplyr::group_by(.data$ID) %>%
        dplyr::count()

      duplicated_n <- duplicated_count_all %>%
        dplyr::filter(.data$n>x) %>%
        dplyr::pull(.data$ID)

      return(duplicated_n)

    }

    over_6 <- repeat_counter(5)
    rel <- rel %>%
      dplyr::filter(!.data$ID1 %in% over_6 & !.data$ID2 %in% over_6)

    over_5 <- repeat_counter(4)
    rel <- rel %>%
      dplyr::filter(!.data$ID1 %in% over_5 & !.data$ID2 %in% over_5)

    over_4 <- repeat_counter(3)
    rel <- rel %>%
      dplyr::filter(!.data$ID1 %in% over_4 & !.data$ID2 %in% over_4)

    over_3 <- repeat_counter(2)
    rel <- rel %>%
      dplyr::filter(!.data$ID1 %in% over_3 & !.data$ID2 %in% over_3)

    over_2 <- repeat_counter(1)
    rel <- rel %>%
      dplyr::filter(!.data$ID1 %in% over_2 & !.data$ID2 %in% over_2)

    exlcusion_over <- c(over_6,over_5,over_4,over_3,over_2)

    remaining_ID <- rel %>%
      dplyr::mutate(excluded_ID = dplyr::case_when(.data$lower_missing==1 ~ .data$ID2,
                                                   is.na(.data$missingness_ID2) ~ .data$ID2,
                                                   .data$lower_missing==2 ~ .data$ID1,
                                                   is.na(.data$missingness_ID1) ~ .data$ID1)) %>%
      dplyr::pull(.data$excluded_ID)
    all_exclusions <- c(remaining_ID,exlcusion_over)

    pheno_new <- pheno %>%
      dplyr::filter(!.data$eid %in% all_exclusions)

    return(pheno_new)

  } else {

    ## now related removed
    cases <- pheno %>%
      dplyr::filter(.data$any_code>0) %>%
      dplyr::pull(.data$eid)

    controls <- pheno %>%
      dplyr::filter(.data$any_code==0) %>%
      dplyr::pull(.data$eid)

    ## need to remove related people for each phenotype
    ## edit to remove exclusion from kinship file
    call_rate_kinship <- call_rate_kinships %>%
      dplyr::filter(.data$ID1 %in% pheno$eid & .data$ID2 %in% pheno$eid)

    ## remove the case with the lowest call rate
    case_case <- call_rate_kinship %>%
      dplyr::filter(.data$ID1 %in% cases, .data$ID2 %in% cases) %>%
      dplyr::mutate(exclusion=dplyr::if_else(.data$lower_missing==1,.data$ID2,.data$ID1)) %>%
      dplyr::select(.data$exclusion) %>%
      dplyr::pull()

    ## edit the callrate_kinship file to remove any excluded cases
    call_rate_kinship <- call_rate_kinship %>%
      dplyr::filter(!.data$ID1 %in% case_case) %>%
      dplyr::filter(!.data$ID2 %in% case_case)

    ## now remove controls, here controls are defined simply as not being exclusions or cases
    case_control <- call_rate_kinship %>%
      dplyr::filter(.data$ID1 %in% cases & .data$ID2 %in% controls) %>%
      dplyr::mutate(exclusion=.data$ID2) %>%
      dplyr::select(.data$exclusion) %>%
      dplyr::pull()

    ## edit the callrate_kinship file to remove any excluded controls
    call_rate_kinship <- call_rate_kinship %>%
      dplyr::filter(!.data$ID1 %in% case_control) %>%
      dplyr::filter(!.data$ID2 %in% case_control)

    control_case <- call_rate_kinship %>%
      dplyr::filter(.data$ID1 %in% controls & .data$ID2 %in% cases) %>%
      dplyr::mutate(exclusion=.data$ID1) %>%
      dplyr::select(.data$exclusion) %>%
      dplyr::pull()

    ## edit the callrate_kinship file to remove any excluded controls
    call_rate_kinship <- call_rate_kinship %>%
      dplyr::filter(!.data$ID1 %in% control_case) %>%
      dplyr::filter(!.data$ID2 %in% control_case)

    ##
    control_control <- call_rate_kinship %>%
      dplyr::filter(.data$ID1 %in% controls & .data$ID2 %in% controls) %>%
      dplyr::mutate(exclusion=dplyr::if_else(.data$lower_missing==1,.data$ID2,.data$ID1)) %>%
      dplyr::select(.data$exclusion) %>%
      dplyr::pull()

    ## combining the vectors
    excluded_cases <- unique(c(case_case))
    excluded_controls <- unique(c(case_control,control_case,control_control))
    removed_ID <- c(excluded_cases,excluded_controls)

    related_removed_final <- pheno %>%
      dplyr::filter(!.data$eid %in% removed_ID)

    return(related_removed_final)
  }
}
#' Creates female specific versions of existing phenotypes
#'
#' @param x Name of sex-specific phenotype
#' @param y Non-sex specific phenotype to be edited
#' @param female_ID vector of female IDs
#' @return Female sex-specific phenotype
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
sex_specific_female <- function(x,y,female_ID) {
  female_pheno <- y %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$eid %in% female_ID) %>%
    dplyr::rename(!!x:=2)
  return(female_pheno)
}
#' Creates male specific versions of existing phenotypes
#'
#' @param x Name of sex-specific phenotype
#' @param y Non-sex specific phenotype to be edited
#' @param male_ID vector of male IDs
#' @return Male sex-specific phenotype
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
sex_specific_male <- function(x,y,male_ID) {
  male_pheno <- y %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$eid %in% male_ID) %>%
    dplyr::rename(!!x:=2)
  return(male_pheno)
}
#' Creates age of onset specific versions of existing phenotypes
#'
#' @param a Phenotype data-frame
#' @param b Non-sex specific phenotype to be edited
#' @param c vector of male IDs
#' @param d Name of age_of_onset phenotype
#' @param DOB_df date of birth data frame
#' @return Age of onset phenotype as list.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
age_of_onset_pheno_creator <- function(d,a,b,c,DOB_df){
  if(is.na(as.numeric(b))){
    rlang::abort(paste0("'lower_limit' must be a numeral"))
  }
  if(is.na(as.numeric(c))){
    rlang::abort(paste0("'upper_limit' must be a numeral"))
  }
  lower_cut_off <- as.numeric(b)
  upper_cut_off <- as.numeric(c)
  age_of_onset_phenotype <- a %>%
    dplyr::left_join(DOB_df, by="eid")%>%
    dplyr::rename(any_code=2) %>%
    dplyr::filter(.data$any_code >=1) %>%
    dplyr::mutate(age_onset=lubridate::interval(.data$DOB, .data$earliest_date) / lubridate::years(1)) %>%
    dplyr::filter(dplyr::between(.data$age_onset,lower_cut_off,upper_cut_off)) %>%
    dplyr::select(.data$eid,!!d:=2, .data$earliest_date)
}
#' Splits phenotype by grouping variable
#'
#' @param x group identifier IDs
#' @param all_phenotypes all phenotypes loaded in
#' @return all phenotypes split by each population
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
spliting_by_population <- function(x,all_phenotypes) {
  per_population <- mapply(per_pheno_ancestry,all_phenotypes,SIMPLIFY = F, USE.NAMES = T,MoreArgs = list(y=x))
  return(per_population)
}
#' Splits phenotype by grouping variable
#'
#' @param x phenotype dataframe
#' @param y Data frame of IDs per group identifier
#' @return phenotype split by population
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
per_pheno_ancestry <- function (x,y) {

  split_pheno <- x %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$eid %in% y)

  return(split_pheno)
}
#' Applies sum_stats per grouping population for each phenotype
#'
#' @param x phenotype dataframe having been split by grouping variable.
#' @param y Grouping variable names.
#' @param PheWAS_manifest manifest file
#' @return summary stats per phenotype per grouping variable.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
summary_all_pops <- function(x,y,PheWAS_manifest) {

  population_name <- y
  pheno_names <- names(x)

  summary_N <- mapply(sum_stats,
                      x,
                      z=pheno_names,
                      MoreArgs = list(y=population_name,PheWAS_manifest=PheWAS_manifest),
                      SIMPLIFY = F, USE.NAMES = T) %>%
    purrr::reduce(dplyr::bind_rows)

  return(summary_N)

}
#' Applies sum_stats per grouping population for each phenotype returning PheWAS_IDs filtered by case number
#'
#' @param x phenotype dataframe having been split by grouping variable.
#' @param quant_min_cases minimum case number required for quantitative traits.
#' @param binary_min_cases minimum case number required for binary traits.
#' @param PheWAS_manifest PheWAS_manifest.
#' @return PheWAS_ID meet filtering conditions for case number.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
selecting_phenotype_cases_N <- function(x,quant_min_cases,binary_min_cases,PheWAS_manifest) {

  PheWAS_ID_type <- PheWAS_manifest %>%
    dplyr::select(.data$PheWAS_ID,.data$analysis)

    filtering_process <- x %>%
      dplyr::rename(cases=2) %>%
      dplyr::mutate(ID_2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female"),stringr::str_remove(.data$PheWAS_ID,"_male|_female"),
                                ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),stringr::str_remove(.data$PheWAS_ID,"_age_of_onset"),.data$PheWAS_ID))) %>%
      dplyr::left_join(PheWAS_ID_type, by=c("ID_2"="PheWAS_ID")) %>%
      dplyr::mutate(analysis=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),"quant",.data$analysis)) %>%
      dplyr::mutate(sufficent_cases=ifelse(.data$analysis=="quant" &
                                             .data$cases >= quant_min_cases,1,ifelse(.data$analysis=="binary" &
                                                                                       .data$cases >= binary_min_cases, 1, 0))) %>%
      dplyr::filter(.data$sufficent_cases>0)

  N_cases_filtered <- x %>%
    dplyr::filter(.data$PheWAS_ID %in% filtering_process$PheWAS_ID) %>%
    dplyr::pull(.data$PheWAS_ID)

  return(N_cases_filtered)
}
#' Removes phenotypes based on case number.
#'
#' @param x all phenotypes as list.
#' @param y PheWAS_IDs to include.
#' @return Phenotypes filtered for case number
#' @keywords internal
reducing_phenotypes <- function(x,y) {

  x <- x[names(x) %in% y == TRUE]
  return(x)
}
#' Removes related individuals from phenotypes per grouping.
#'
#' @param x all phenotypes as list.
#' @param call_rate_kinships call_rate_kinship file
#' @param PheWAS_manifest manifest file
#' @return Phenotypes filtered for relatedness.
#' @keywords internal
relate_remove_all_pops <- function(x,call_rate_kinships,PheWAS_manifest) {

  PheWAS_ID <- names(x)
  relate_remove_step <- mapply(relate_remove_fn,x,
                               PheWAS_ID,
                               MoreArgs = list(call_rate_kinships=call_rate_kinships,PheWAS_manifest=PheWAS_manifest),
                               SIMPLIFY = F,USE.NAMES = T)
  return(relate_remove_step)
}

#' Allows main conversion to data frame of phenotypes
#'
#' @param x All phenotypes as list objects
#' @param IVNT logical value if using IVNT transformations
#' @param age_of_onset_phenotypes data frame of age_of_onset phenotypes or a NULL value
#' @param PheWAS_manifest manifest document
#' @param age_of_onset_all if True creates all age of onset phenotypes
#' @param all_phenotypes all available phenotypes
#' @return Phenotypes as dataframe ready for analysis
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
join_cols_pheno <- function (x,IVNT,age_of_onset_phenotypes,PheWAS_manifest,age_of_onset_all,all_phenotypes) {
  per_list <- x
  PheWAS_ID <- names(per_list)
  final_phenotypes <- mapply(converted_for_association,
                             per_list,
                             PheWAS_ID,
                             MoreArgs = list(PheWAS_manifest=PheWAS_manifest,
                                             PheWAS_ID_list=PheWAS_ID),
                             SIMPLIFY = F,USE.NAMES = T)

  if (IVNT) {

    PheWAS_IDs <- names(final_phenotypes)
    transformed_phenoptypes_final <- mapply(IVNT_transformation,
                                            final_phenotypes,
                                            PheWAS_IDs,
                                            MoreArgs = list(PheWAS_manifest=PheWAS_manifest,
                                                            age_of_onset_phenotypes=age_of_onset_phenotypes,
                                                            age_of_onset_all=age_of_onset_all,
                                                            all_phenotypes=all_phenotypes),
                                            SIMPLIFY = F,USE.NAMES = T) %>%
      purrr::reduce(dplyr::full_join)

    return(transformed_phenoptypes_final)

  } else {

    final_phenotypes_tabled <- final_phenotypes %>%
      purrr::reduce(dplyr::full_join)

    return(final_phenotypes_tabled)

  }

}
#' Converts phenotypes to form suitable for association analysis retains list format
#'
#' @param x Phenotypes as list object
#' @param y PheWAS_ID of phenotype
#' @param PheWAS_manifest the PheWAS manifest
#' @param PheWAS_ID_list all PheWAS ID's in analysis
#' @return Phenotypes as list edited to for combination into dataframe for association analysis
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
converted_for_association <- function (x,y,PheWAS_manifest,PheWAS_ID_list) {

  if(stringr::str_detect(y,"age_of_onset")){
    analysis_type <- "quant"
    age_col <- "age"
  } else {
  name <- stringr::str_remove(y,"_male|_female")
  analysis_type <- PheWAS_manifest %>%
    dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
    dplyr::select(.data$analysis) %>%
    dplyr::pull()

  age_col <- PheWAS_manifest %>%
    dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
    dplyr::select(.data$age_col) %>%
    dplyr::pull()
  }

  if (analysis_type == "binary") {
    phenotype <- x %>%
      dplyr::rename(any_code=2) %>%
      dplyr::mutate(!!y:=ifelse(is.na(.data$any_code),NA,ifelse(.data$any_code>0,1,0))) %>%
      dplyr::select(.data$eid,{{y}}) %>%
      tidyr::drop_na()
    return(phenotype)

  } else if (analysis_type == "quant"){
    if(stringr::str_detect(y,"age_of_onset")){
      desired_col_names <- c("eid",y)
    }
    if(stringr::str_detect(y,"_male|_female") & stringr::str_remove(y,"_male|_female") %in% PheWAS_ID_list) {
      desired_col_names <- c("eid",y)
    } else {
      if(age_col=="named") {
        desired_col_names <- c("eid",y,paste0(name,"_age"))
      } else {
        desired_col_names <- c("eid",y)
      }
    }
    phenotype <- x %>%
      dplyr::rename(!!y:=2) %>%
      dplyr::select(tidyselect::any_of(desired_col_names)) %>%
  tidyr::drop_na()

return(phenotype)
  }
}
#' Applies transformations to phenotypes as directed by manifest file.
#'
#' @param x Phenotypes as list object
#' @param y PheWAS_ID of phenotype
#' @param PheWAS_manifest the PheWAS manifest
#' @param age_of_onset_phenotypes either NULL or location of age_of_onset phenotypes file
#' @param age_of_onset_all If T creates all age_of_onset phenotypes
#' @param all_phenotypes all available phenotypes
#' @return Phenotypes as list transformed as directed ready to convert into a dataframe
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
IVNT_transformation <- function(x,y,PheWAS_manifest,age_of_onset_phenotypes,age_of_onset_all,all_phenotypes) {
  if(any(!is.null(age_of_onset_phenotypes) | age_of_onset_all)){
    if(!is.null(age_of_onset_phenotypes)){
  age_of_onset_df <- data.table::fread(age_of_onset_phenotypes)
  }
  else{
  age_of_onset_df <- PheWAS_manifest %>%
    dplyr::filter(.data$included_in_analysis == 1) %>%
    dplyr::filter(.data$analysis == "binary") %>%
    dplyr::select(.data$PheWAS_ID) %>%
    dplyr::mutate(lower_limit=1,upper_limit=120,transformation="IVNT") %>%
    dplyr::filter(.data$PheWAS_ID %in% names(all_phenotypes))
  }
}
  if(stringr::str_detect(y,"age_of_onset")){
    name <- stringr::str_remove(y,"_age_of_onset")
    type <- age_of_onset_df %>%
      dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
      dplyr::pull(.data$transformation)
    age_col <- "age"
  } else {
  name <- stringr::str_remove(y,"_male|_female")

  age_col <- PheWAS_manifest %>%
    dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
    dplyr::select(.data$age_col) %>%
    dplyr::pull()

  type <- PheWAS_manifest %>%
    dplyr::filter(.data$PheWAS_ID=={{name}}) %>%
    dplyr::select(.data$transformation) %>%
    dplyr::pull()
}
  if (type == "IVNT") {

    pheno <- x %>%
      dplyr::rename(any_code=2)

    pheno <- data.table::as.data.table(pheno)

    pheno_col <- pheno %>%
      dplyr::pull(.data$any_code)

    message(y)

    IVNT_calc <- RNOmni::RankNorm(pheno_col)

    IVNT_col <- replace(pheno[,2],values = IVNT_calc) %>%
      dplyr::rename(!!y:=1)

    if(age_col=="named"){
      desired_col_names <- c("eid",y,paste0(name,"_age"))
    } else {
      desired_col_names <- c("eid",y)
    }
    pheno_new <- pheno %>%
      dplyr::bind_cols(IVNT_col) %>%
      dplyr::select(tidyselect::any_of(desired_col_names))

    return(pheno_new)
  } else {
      return(x)
  } 
}
#' Saves list objects as dataframes per grouping variable
#'
#' @param a Phenotypes as list object
#' @param b grouping name variable
#' @param save_name_core save_name given as argument
#' @param save_folder folder for the save_name given as an argument.
#' @return A saved table of phenotypes per grouping name variable.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
saving_files <- function(a,b,save_name_core,save_folder) {

  phenos_table <- a
  if(stringr::str_detect(save_name_core,".gz")){
    name <- paste0(b,"_",save_name_core)
  } else{
    name <- paste0(b,"_",save_name_core,".gz")
  }
  data.table::fwrite(phenos_table,paste0(save_folder,"/",name),sep = "\t",na = NA)
}

#' Filters and prepares phenotyeps for association analysis.
#'
#' @param phenotype_filtered_save_name Full file path for save location of filtered phenotypes. Saves and dataframe per-group (if provided) that is the input to the association testing. Name of final save is /path/to/file/(group)_savename.
#' @param phenotype_folder Full path of the folder containing the phenotype data created in previous steps.
#' @param phenotype_files Comma separated full file paths of phenotype data created in previous steps.
#' @param groupings Full path of the file containing group information used for stratifying the sample. If this argument is specified, the phenotype filtering will be performed within the groups. Otherwise, the filtering will be performed in all individuals.
#' @param quantitative_Case_N Number that represents the minimum number of cases for quantitative phenotype inclusion default=100.
#' @param binary_Case_N Number that represents the minimum number of cases for binary phenotype inclusion default=50.
#' @param relate_remove Specify whether related individuals are additionally excluded (TRUE) or not (FALSE). Default is FALSE.
#' @param kinship_file Full path to the file containing kinship coefficient estimates.
#' @param sex_split_all Specify whether all phenotypes should be made into sex specific phenotypes (TRUE) or not (FALSE). Only one of sex_split_all and sex_split_phenotypes, can be selected all sex specific phenotypes be labelled (PheWAS_ID)_male or (PheWAS_ID)_female.. Default to FALSE.
#' @param sex_split_phenotypes Full file path to file containing single column of PheWAS_IDs to make sex specific phenotypes, will derive and add a male and female version of each phenotype. All sex specific phenotypes be labelled (PheWAS_ID)_male or (PheWAS_ID)_female.
#' @param sex_info Full path of the combined_sex file (generated by the 02_data_preparation.R script).
#' @param male Value corresponding to males in the combined_sex file must be a number. default=1.
#' @param female Value corresponding to females in the combined_sex file must be a number. default=0.
#' @param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#' @param age_of_onset_phenotypes Full file path to file containing three columns PheWAS_ID,lower_limit,upper_limit. PheWAS_ID are the phenotypes to create the age_of_onset phenotypes for, lower_limit is the lower age boundary to filter read as <=, upper_limit is the upper age boundary acceptable read as >=. All age of onset phenotypes to be labelled (PheWAS_ID)_age_of_onset.
#' @param age_of_onset_all Specify if wanting to create an age of onset version for all suitable phenotypes (TRUE) will default lower limit to 30 and upper limit to 120 and request IVNT transformation for each created phenotype. Default is FALSE.
#' @param DOB_file Full file path to file containing DOB information.
#' @param IVNT Specify whether quantitative phenotypes, defined in the PheWAS manifest, should undergo a rank-based inverse normal transformation (TRUE) or not (FALSE). Default is FALSE.
#' @param save_RDS Full file path to save location to save and intermediate RDS file that has been filtered by case_N and relatedness (if chosen) but not converted to dataframe for analysis.saves an intermediate RDS file that has been filtered by case_N and relatedness (if chosen) but not converted to dataframe for analysis.
#' @param stats_save Full file path to save file for stats, saves summary stats for phenotypes for number of cases and controls in each grouping variable.File name will be appended with _N_filtered.csv or _relate_remove.csv, depending on if relate remove was selected.
#' @return A data frame per grouping variable of prepared and filtered phenotypes used for association analysis.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

phenotype_preparation <- function(phenotype_filtered_save_name,
                                  phenotype_folder,
                                  phenotype_files,
                                  groupings,
                                  quantitative_Case_N,
                                  binary_Case_N,
                                  relate_remove,
                                  kinship_file,
                                  sex_split_all,
                                  sex_split_phenotypes,
                                  sex_info,
                                  male,
                                  female,
                                  PheWAS_manifest_overide,
                                  age_of_onset_phenotypes,
                                  age_of_onset_all,
                                  DOB_file,
                                  IVNT,
                                  save_RDS,
                                  stats_save){

# load in
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

if(!is.null(phenotype_folder)){
  phenotype_files_vector <- list.files(phenotype_folder,pattern = ".RDS",full.names = T)
} else if(!is.null(phenotype_files)){
  phenotype_files_vector <- c(unlist(strsplit(phenotype_files,",")))
}

if(!dir.exists(dirname(phenotype_filtered_save_name))){
  dir.create(dirname(phenotype_filtered_save_name))}
save_name_core <- basename(phenotype_filtered_save_name)
save_folder <- dirname(phenotype_filtered_save_name)

if(is.na(as.numeric(quantitative_Case_N))){
  rlang::abort(paste0("'quantitative_Case_N' must be a numeral"))
}
quant_min_cases <- as.numeric(quantitative_Case_N)

if(is.na(as.numeric(binary_Case_N))){
  rlang::abort(paste0("'binary_Case_N' must be a numeral"))
}
binary_min_cases <- as.numeric(binary_Case_N)

# load optional arguments
if(!is.null(groupings)) {
  if(!file.exists(groupings)){
    rlang::abort(paste0("'groupings' must be a file"))
  }
  groupings_df <- data.table::fread(groupings)
  groupings_colname <- c("eid","group")

  if(!all(tibble::has_name(groupings_df,groupings_colname))){
    warning(paste0("'groupings' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(groupings_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(groupings_df), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(groupings_df),groupings_colname), collapse=","))

  }
}

# read in all the files in the phenotypes folder and combine them.
phenotypes <- phenotype_files_vector %>%
  purrr::map(readRDS) %>%
  purrr::reduce(c)

# first edit the phenotypes to those included in the phenotype manifest for inclusion
included_phenotypes <- PheWAS_manifest %>%
  dplyr::filter(.data$included_in_analysis == 1)
all_phenotypes <- phenotypes[names(phenotypes) %in% included_phenotypes$PheWAS_ID == TRUE]

# optional addition of sex_split phenotypes

if(any(!is.null(sex_split_phenotypes) | sex_split_all)){

  pre_sex_split <- PheWAS_manifest %>%
    dplyr::filter(.data$sex=="Male"|.data$sex=="Female")

  all_phenotypes_minus_existing_sexpheno<-all_phenotypes[names(all_phenotypes) %in% pre_sex_split$PheWAS_ID == F]

  if(!is.null(sex_split_phenotypes)){

  if(!file.exists(sex_split_phenotypes)){
    rlang::abort(paste0("'sex_split_phenotypes' must be a file"))
  }

  sex_phenos <- data.table::fread(sex_split_phenotypes)
  sex_phenos_colname <- c("PheWAS_ID")

  if(!all(tibble::has_name(sex_phenos,sex_phenos_colname))){
    warning(paste0("'sex_split_phenotypes' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(sex_phenos_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(sex_phenos), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(sex_phenos),sex_phenos_colname), collapse=","))

  }
  } else {
    sex_phenos <- data.frame(PheWAS_ID=names(all_phenotypes_minus_existing_sexpheno))
  }
  if(is.na(as.numeric(male))){
    rlang::abort(paste0("'male' must be a numeral"))
  }
  male_N <- as.numeric(male)
  if(is.na(as.numeric(female))){
    rlang::abort(paste0("'female' must be a numeral"))
  }
  female_N <- as.numeric(female)

  if(!file.exists(sex_info)){
    rlang::abort(paste0("'sex_info' must be a file"))
  }
  sex_info_df <- data.table::fread(sex_info)
  sex_info_colname <- c("eid","sex")

  if(!all(tibble::has_name(sex_info_df,sex_info_colname))){
    warning(paste0("'sex_info' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(sex_info_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(sex_info_df), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(sex_info_df),sex_info_colname), collapse=","))

  }

  male_ID <-  sex_info_df%>%
    dplyr::filter(.data$sex==male_N) %>%
    dplyr::pull(.data$eid)
  female_ID <- sex_info_df %>%
    dplyr::filter(.data$sex==female_N) %>%
    dplyr::pull(.data$eid)

  all_phenotypes_sex_edit <- all_phenotypes_minus_existing_sexpheno[names(all_phenotypes_minus_existing_sexpheno) %in% sex_phenos$PheWAS_ID == T]

  male_names <- paste0(names(all_phenotypes_sex_edit),"_male")
  female_names <- paste0(names(all_phenotypes_sex_edit),"_female")

  female_phenotypes <- mapply(sex_specific_female,
                              female_names,
                              all_phenotypes_sex_edit,
                              MoreArgs = list(female_ID=female_ID),
                              SIMPLIFY = F, USE.NAMES = T)
  male_phenotypes <- mapply(sex_specific_male,
                            male_names,
                            all_phenotypes_sex_edit,
                            MoreArgs = list(male_ID=male_ID),
                            SIMPLIFY = F, USE.NAMES = T)

  all_phenotypes <- c(all_phenotypes,female_phenotypes,male_phenotypes)

}

# age of onset creation
if(any(!is.null(age_of_onset_phenotypes) | age_of_onset_all)){

  if(!is.null(age_of_onset_phenotypes)){
  if(!file.exists(age_of_onset_phenotypes)){
    rlang::abort(paste0("'age_of_onset_phenotypes' must be a file"))
  }
  age_of_onset <- data.table::fread(age_of_onset_phenotypes)
  age_of_onset_colname <- c("PheWAS_ID","lower_limit","upper_limit","transformation")

  if(!all(tibble::has_name(age_of_onset,age_of_onset_colname))){
    warning(paste0("'age_of_onset_phenotypes' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(age_of_onset_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(age_of_onset), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(age_of_onset),age_of_onset_colname), collapse=","))

  }
  } else {
    age_of_onset <- PheWAS_manifest %>%
      dplyr::filter(.data$included_in_analysis == 1) %>%
      dplyr::filter(.data$analysis == "binary") %>%
      dplyr::select(.data$PheWAS_ID) %>%
      dplyr::mutate(lower_limit=1,upper_limit=120,transformation="IVNT") %>%
      dplyr::filter(.data$PheWAS_ID %in% names(all_phenotypes))

}
  if(!file.exists(DOB_file)){
    rlang::abort(paste0("'DOB_file' must be a file"))
  }
  DOB_df <- data.table::fread(DOB_file)
  DOB_df_colname <- c("eid","DOB")

  if(!all(tibble::has_name(DOB_df,DOB_df_colname))){
    warning(paste0("'DOB_file' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(DOB_df_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(DOB_df), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(DOB_df),DOB_df_colname), collapse=","))

  }

  all_phenotypes_AOO <- all_phenotypes[names(all_phenotypes) %in% age_of_onset$PheWAS_ID == T]

  pheno_names <- paste0(names(all_phenotypes_AOO),"_age_of_onset")

  AOO_phenotypes <- mapply(age_of_onset_pheno_creator,
                           pheno_names,
                           all_phenotypes_AOO,
                           age_of_onset$lower_limit,
                           age_of_onset$upper_limit,
                           MoreArgs = list(DOB_df=DOB_df),
                           SIMPLIFY = F, USE.NAMES = T)

  all_phenotypes <- c(all_phenotypes,AOO_phenotypes)
}
# now split the phenotypes into ancestry or any other division
if(!is.null(groupings)){

  groups <- purrr::discard(unique(groupings_df$group),is.na)
  group_ID <- mapply(group_id_creater,
                     groups,
                     MoreArgs = list(groupings_df=groupings_df),
                     USE.NAMES = T,SIMPLIFY = F)

  per_population_phenotypes <- mapply(spliting_by_population,
                                      group_ID,
                                      MoreArgs = list(all_phenotypes=all_phenotypes),
                                      SIMPLIFY = F,USE.NAMES = T)

} else {
  groups <- c("all_population")
  per_population_ungroup <- all_phenotypes
  per_population_phenotypes <- list(all_pop=per_population_ungroup)
}

# now run summary stats on those phenotypes using whatever groupings selected here we use ancestry
per_population_summary <- mapply(summary_all_pops,
                                 per_population_phenotypes,
                                 groups,
                                 MoreArgs = list(PheWAS_manifest=PheWAS_manifest),
                                 SIMPLIFY = F,USE.NAMES = T)
if(!is.null(stats_save)){
  if(!dir.exists(dirname(stats_save))){
    dir.create(dirname(stats_save))}
  stats_save_df <- per_population_summary %>%
    purrr::reduce(dplyr::full_join) %>%
    dplyr::mutate(ID_2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female"),stringr::str_remove(.data$PheWAS_ID,"_male|_female"),
                              ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),stringr::str_remove(.data$PheWAS_ID,"_age_of_onset"),.data$PheWAS_ID))) %>%
    dplyr::left_join(PheWAS_manifest %>% dplyr::select(.data$PheWAS_ID,.data$phenotype),by=c("ID_2"="PheWAS_ID")) %>%
    dplyr::mutate(phenotype=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),paste0("Age of onset of ",.data$phenotype),
                                   ifelse(stringr::str_detect(.data$PheWAS_ID,"_male"),paste0("Male stratified variant of ",.data$phenotype),
                                          ifelse(stringr::str_detect(.data$PheWAS_ID,"_female"),paste0("Female stratified variant of ",.data$phenotype),
                                                 .data$phenotype)))) %>%
    dplyr::relocate(.data$phenotype,.after = .data$PheWAS_ID) %>%
    dplyr::select(-.data$ID_2)

  save_name_stats <- paste0(stats_save,"_N_filtered.csv")

  data.table::fwrite(stats_save_df,save_name_stats)

} else {

}

# now exclude those with <N cases <N cases quantitative.
trimmed_PheWAS_ID <- mapply(selecting_phenotype_cases_N,
                            per_population_summary,
                            MoreArgs = list(quant_min_cases=quant_min_cases,
                                            binary_min_cases=binary_min_cases,
                                            PheWAS_manifest=PheWAS_manifest),
                            SIMPLIFY = F,USE.NAMES = T)

epi_case_trimmed <- mapply(reducing_phenotypes,per_population_phenotypes,trimmed_PheWAS_ID,SIMPLIFY = F)

# now to remove related individuals within each population.
if(relate_remove) {
  if(!file.exists(kinship_file)){
    rlang::abort(paste0("'kinship_file' must be a file"))
  }
  call_rate_kinships <- data.table::fread(kinship_file)
  call_rate_kinships_colname <- c("ID1","ID2","Kinship","missingness_ID1","missingness_ID2","lower_missing")

  if(!all(tibble::has_name(call_rate_kinships,call_rate_kinships_colname))){
    warning(paste0("'kinship_file' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(call_rate_kinships_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(call_rate_kinships), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(call_rate_kinships),call_rate_kinships_colname), collapse=","))

  }
  #run functions to remove related individuals from each phenotype
  relate_remove_phenotypes <- mapply(relate_remove_all_pops,
                                     epi_case_trimmed,
                                     MoreArgs = list(call_rate_kinships=call_rate_kinships,
                                                     PheWAS_manifest=PheWAS_manifest),
                                     SIMPLIFY = F,USE.NAMES = T)

  relate_remove_phenotypes <- Filter(length,relate_remove_phenotypes)

  RR_population_summary <- mapply(summary_all_pops,
                                  relate_remove_phenotypes,
                                  names(relate_remove_phenotypes),
                                  MoreArgs = list(PheWAS_manifest=PheWAS_manifest),
                                  SIMPLIFY = F,USE.NAMES = T)
  if(!is.null(stats_save)){
    stats_save_RR <- RR_population_summary %>%
      purrr::reduce(dplyr::full_join) %>%
      dplyr::mutate(ID_2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female"),stringr::str_remove(.data$PheWAS_ID,"_male|_female"),
                                ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),stringr::str_remove(.data$PheWAS_ID,"_age_of_onset"),.data$PheWAS_ID))) %>%
      dplyr::left_join(PheWAS_manifest %>% dplyr::select(.data$PheWAS_ID,.data$phenotype),by=c("ID_2"="PheWAS_ID")) %>%
      dplyr::mutate(phenotype=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),paste0("Age of onset of ",.data$phenotype),
                                     ifelse(stringr::str_detect(.data$PheWAS_ID,"_male"),paste0("Male stratified variant of ",.data$phenotype),
                                            ifelse(stringr::str_detect(.data$PheWAS_ID,"_female"),paste0("Female stratified variant of ",.data$phenotype),
                                                   .data$phenotype)))) %>%
      dplyr::relocate(.data$phenotype,.after = .data$PheWAS_ID) %>%
      dplyr::select(-.data$ID_2)

    save_name_stats_RR <- paste0(stats_save,"_relate_remove.csv")

    data.table::fwrite(stats_save_RR,save_name_stats_RR)

  } else {

  }
  RR_trimmed_PheWAS_ID <- mapply(selecting_phenotype_cases_N,
                                 RR_population_summary,
                                 MoreArgs = list(quant_min_cases=quant_min_cases,
                                                 binary_min_cases=binary_min_cases,
                                                 PheWAS_manifest=PheWAS_manifest),
                                 SIMPLIFY = F,USE.NAMES = T)

  RR_case_trimmed <- mapply(reducing_phenotypes,
                            relate_remove_phenotypes,
                            RR_trimmed_PheWAS_ID,
                            SIMPLIFY = F)

  pheno_list_final <- RR_case_trimmed

} else {
  pheno_list_final <- epi_case_trimmed
}

if(!is.null(save_RDS)){
  if(!dir.exists(dirname(save_RDS))){
    dir.create(dirname(save_RDS))
  }
  saveRDS(pheno_list_final,save_RDS)
}

final_phenotypes_columns <- mapply(join_cols_pheno,
                                   pheno_list_final,
                                   MoreArgs = list(PheWAS_manifest=PheWAS_manifest,
                                                   IVNT=IVNT,
                                                   age_of_onset_phenotypes=age_of_onset_phenotypes,
                                                   age_of_onset_all=age_of_onset_all,
                                                   all_phenotypes=all_phenotypes),
                                   SIMPLIFY = F,USE.NAMES = T)

mapply(saving_files,
       final_phenotypes_columns,
       names(final_phenotypes_columns),
       MoreArgs = list(save_name_core=save_name_core,
                       save_folder=save_folder),
       SIMPLIFY = F)
}

