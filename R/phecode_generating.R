#' Creates phecode phenotypes for available ICD9 and ICD10 data
#'
#' @param a PheWAS_ID
#' @param b phecode
#' @param c 1st lower phecode range
#' @param d 1st upper phecode range
#' @param e 2nd lower phecode range
#' @param f 2nd upper phecode range
#' @param g 3rd lower phecode range
#' @param h 3rd upper phecode range
#' @param i if sex specific or not
#' @param j data frame of mapped phecodes to ICD10/9 data
#' @param k data frame for control exclusions
#' @param l data frame fro combined sex
#' @return A dataframe containing phecode status for each participant
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom rlang .data

phecoding <- function(a,b,c,d,e,f,g,h,i,j,k,l) {
  # j is sex, a small number of phecodes are sex specific and so the function needs to account for this.
  if (is.na(i)) {
    # searching for data
    phecodes <- j %>%
      dplyr::filter(.data$phecode==b)
    # making cases
    cases_epi_date <- phecodes %>%
      dplyr::group_by(.data$eid) %>%
      dplyr::summarise(earliest_date=min(.data$date))
    cases_epi <- phecodes %>%
      dplyr::count(.data$eid,.data$phecode) %>%
      dplyr::left_join(cases_epi_date) %>%
      dplyr::mutate(!!a := .data$n) %>%
      dplyr::select(.data$eid,{{a}},.data$earliest_date)
    case_id_epi <- cases_epi %>%
      dplyr::select(1) %>%
      dplyr::pull()
    #excluded are those in the exclusion ranges so these cannot be controls.
    excluded_epi_df <- j %>%
      dplyr::filter(!.data$eid %in% case_id_epi) %>%
      dplyr::filter(dplyr::between(.data$phecode,c,d) | dplyr::between(.data$phecode,e,f) | dplyr::between(.data$phecode,g,h)) %>%
      dplyr::mutate(!!a := NA,
                    earliest_date=NA,
                    earliest_date=lubridate::ymd(.data$earliest_date)) %>%
      dplyr::distinct(.data$eid, .keep_all = T) %>%
      dplyr::select("eid",{{a}},"earliest_date")
    # vector to use later
    excluded_epi <- excluded_epi_df %>%
      dplyr::select(.data$eid) %>%
      dplyr::pull()
    #controls are then dplyr::selected from the remaining population sex label has virtually everyone in UK Biobank and so is the population
    #sample for the epidemiological definition
    controls_epi <- l %>%
      dplyr::filter(!.data$eid %in% case_id_epi) %>%
      dplyr::filter(!.data$eid %in% excluded_epi) %>%
      dplyr::filter(!.data$eid %in% k) %>%
      dplyr::mutate(!!a := 0,
                    earliest_date =NA,
                    earliest_date=lubridate::ymd(.data$earliest_date)) %>%
      dplyr::select(.data$eid,{{a}},.data$earliest_date)
    # final output
    phecode_status <- cases_epi %>%
      dplyr::mutate(earliest_date=as.Date(.data$earliest_date)) %>%
      dplyr::bind_rows(controls_epi, excluded_epi_df)
    return(phecode_status)
  } else {
    #exactly as above but dplyr::filtering for sex
    phecodes <- j %>%
      dplyr::filter(.data$sex==i) %>%
      dplyr::filter(.data$phecode==b)
    # making cases
    cases_epi_date <- phecodes %>%
      dplyr::group_by(.data$eid) %>%
      dplyr::summarise(earliest_date=min(.data$date))
    cases_epi <- phecodes %>%
      dplyr::count(.data$eid,.data$phecode) %>%
      dplyr::left_join(cases_epi_date) %>%
      dplyr::mutate(!!a := .data$n) %>%
      dplyr::select(.data$eid,{{a}},.data$earliest_date)
    # vector used later
    case_id_epi <- cases_epi %>%
      dplyr::select(1) %>%
      dplyr::pull()
    # now exclusions
    excluded_epi_df <- j %>%
      dplyr::filter(.data$sex==i) %>%
      dplyr::filter(dplyr::between(.data$phecode,c,d) | dplyr::between(.data$phecode,e,f) | dplyr::between(.data$phecode,g,h)) %>%
      dplyr::mutate(!!a := NA,
                    earliest_date=NA,
                    earliest_date=lubridate::ymd(.data$earliest_date)) %>%
      dplyr::distinct(.data$eid, .keep_all = T) %>%
      dplyr::select(.data$eid,{{a}},.data$earliest_date)
    # vector used later
    excluded_epi <- excluded_epi_df %>%
      dplyr::select(.data$eid) %>%
      dplyr::pull()
    #controls are then dplyr::selected from the remaining population sex label has virtually everyone in UK Biobank and so is the population
    #sample for the epidemiological definition
    controls_epi <- l %>%
      dplyr::filter(.data$sex==i) %>%
      dplyr::filter(!.data$eid %in% case_id_epi) %>%
      dplyr::filter(!.data$eid %in% excluded_epi) %>%
      dplyr::filter(!.data$eid %in% k) %>%
      dplyr::mutate(!!a := 0,
                    earliest_date =NA,
                    earliest_date=lubridate::ymd(.data$earliest_date)) %>%
      dplyr::select(.data$eid,{{a}},.data$earliest_date)
    # final output
    phecode_status <- cases_epi %>%
      dplyr::mutate(earliest_date=as.Date(.data$earliest_date)) %>%
      dplyr::bind_rows(controls_epi, excluded_epi_df)
    return(phecode_status)
  }
}

#' Creates exclusion ranges which are then used for composite phenotypes
#'
#' @param a Name of range exclusion file
#' @param b 1st lower phecode range
#' @param c 1st upper phecode range
#' @param d 2nd lower phecode range
#' @param e 2nd upper phecode range
#' @param f 3rd lower phecode range
#' @param g 3rd upper phecode range
#' @param h phecodes_mapped data frame
#' @return A dataframe containing ids of participants to be excluded as controls for given set of phenotypes
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
exclusion_generator <- function(a,b,c,d,e,f,g,h) {
  #use the largest data file which includes the mortality and cancer data.
  exclusions <- h %>%
    dplyr::filter(dplyr::between(.data$phecode,b,c) | dplyr::between(.data$phecode,d,e) | dplyr::between(.data$phecode,f,g)) %>%
    dplyr::mutate(!!a := 1) %>%
    dplyr::distinct(.data$eid)
  return(exclusions)
}

#' Creates phecode phenotypes and saves as a list alongside range_exclusions
#'
#' @param health_data Full path of the file produced by the previous step (phenotype_data.R).
#' @param sex_info Full file path of the combined_sex file a file containing participant ID and sex information (0=female, 1=male).
#' @param no_range_ID Choose whether to save lists of participant identifiers to be excluded from control definitions (FALSE) or not (TRUE). Defaults to FALSE.
#' @param no_phecodes Choose whether to save the Phecode-generated phenotype data (FALSE) or not (TRUE). Defaults to FALSE.
#' @param control_exclusions Full file path of the optional control exclusions file.
#' @param N_cores Number of cores requested if parallel computing is desired. Defaults to single core computing.
#' @param phecode_save_file Full path of the save file for the Phecode phenotypes R data object.
#' @param range_ID_save_file Full path of the folder used to store the control exclusions R data object.
#' @param ICD10 Comma-separated string representing the column names used in health_data for ICD-10 values. If there are no ICD-10 values, use NA or any text not used as a source in health data. Defaults to "ICD10".
#' @param ICD9 Comma-separated string representing the column names used in health_data for ICD-9 values. If there are no ICD-9 values, use NA or any text not used as a source in health data. Defaults to "ICD9".
#' @param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#' @return A rds file containing phecode phenotypes.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data

generating_phecodes <- function(health_data,
                                sex_info,
                                no_range_ID,
                                no_phecodes,
                                control_exclusions,
                                N_cores,
                                phecode_save_file,
                                range_ID_save_file,
                                ICD10,
                                ICD9,
                                PheWAS_manifest_overide) {

  # Load in and define variables -----------------------------------------------------------------
  if(!is.null(N_cores)) {
    if(is.na(as.numeric(N_cores))){
      rlang::abort(paste0("'N_cores' must be a numeral"))
    }
    N_cores <- as.numeric(N_cores)
    if(N_cores > parallel::detectCores()){
      rlang::abort(paste0("The number of cores detected is ",parallel::detectCores()," which is lower than the requested 'N_core' of ",N_cores," please lower the requested N_cores to be <= ",parallel::detectCores()))
    }
  } else {
    N_cores <- NA
  }

  phecode_save_location <- phecode_save_file
  new_folder <- stringr::str_remove(phecode_save_file,basename(phecode_save_file))
  if(!dir.exists(new_folder)){
    dir.create(new_folder)
  }

  range_ID_save_location <- range_ID_save_file
  new_folder_range <- stringr::str_remove(range_ID_save_file,basename(range_ID_save_file))
  if(!dir.exists(new_folder_range)){
    dir.create(new_folder_range)}

  # PheWAS manifest
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

  # read in the phecode mapping files for ICD10 and ICD9
  phecode_map_ICD10 <- data.table::fread(system.file("extdata","phecode_map_rollup_ICD10.gz", package = "DeepPheWAS"))
  phecode_map_ICD9 <- data.table::fread(system.file("extdata","phecode_map_rollup_ICD9.gz", package = "DeepPheWAS"))

  # phecode definitions for mapping to phecodes
  phecode_definition <- data.table::fread(system.file("extdata","phecode_definitions.gz", package = "DeepPheWAS"))

  # optional file for control exclusions due to lack of suitable follow-up
  if(is.null(control_exclusions)) {
    control_exclusions <- data.frame(list("eid"="")) %>%
      dplyr::pull(eid)
  } else {
    if(!file.exists(control_exclusions)){
      rlang::abort(paste0("'control_exclusions' must be a file"))
    }
    control_exclusions <- data.table::fread(control_exclusions) %>%
      dplyr::select(1) %>%
      dplyr::pull()
  }

  # cases and controls can only be defined where recorded sex is available
  if(!file.exists(sex_info)){
    rlang::abort(paste0("'sex_info' must be a file"))
  }
  combined_sex <- data.table::fread(sex_info)
  combined_sex_colname <- c("eid","sex")

  if(!all(tibble::has_name(combined_sex,combined_sex_colname))){
    warning(paste0("'sex_info' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(combined_sex_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(combined_sex), collapse=","),
                  paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(combined_sex),combined_sex_colname), collapse=","))

  }

  # health data
  if(!file.exists(health_data)){
    rlang::abort(paste0("'health_data' must be a file"))
  }
  health_data <- data.table::fread(health_data)
  health_data_colname <- c("eid","code","date","source")

  if(!all(tibble::has_name(health_data,health_data_colname))){
    warning(paste0("'health_data' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(health_data_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(health_data), collapse=","),
                  paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(health_data),health_data_colname), collapse=","))

  }

  # Map to hospital data to phecodes -----------------------------------------------------
  # source names for ICD10 records and ICD9 records respectively
  if(is.null(ICD10)) {
    ICD10 <- c("ICD10_1","ICD10_2","ICD10_3","MD","cancer")
  } else {
    ICD10 <- c(unlist(strsplit(ICD10,",")))
  }
  if(is.null(ICD9)) {
    ICD9 <- c("ICD9_1","ICD9_2","ICD9_3")
  } else {
    ICD9 <- c(unlist(strsplit(ICD9,",")))
  }
  # creating an ICD9 and ICD10 specific records for searching
  all_ICD10 <- health_data %>%
    dplyr::filter(source %in% ICD10)
  all_ICD9 <- health_data %>%
    dplyr::filter(source %in% ICD9)
  ## map to the phecodes separately for the different data sources ICD10 and 9 (phecodes have different maps)
  ICD10_mapped <- dplyr::left_join(all_ICD10,phecode_map_ICD10, by=c("code"="ICD10")) %>%
    dplyr::select(.data$eid,phecode=.data$PHECODE,.data$date) %>%
    tidyr::drop_na()
  ICD9_mapped <- dplyr::left_join(all_ICD9,phecode_map_ICD9, by=c("code"="icd9")) %>%
    dplyr::select(.data$eid,.data$phecode,.data$date) %>%
    tidyr::drop_na()
  ## then join the ICD10 and ICD9 data together and add in sex information via combined sex file
  phecodes_mapped <- dplyr::bind_rows(ICD10_mapped,ICD9_mapped) %>%
    dplyr::left_join(combined_sex)

  # filter the phecode definition by available phecodes in the data
  manifest_phecodes <- PheWAS_manifest %>%
    dplyr::filter(.data$category=="Phecode") %>%
    dplyr::mutate(phecode=stringr::str_remove(.data$PheWAS_ID,"P"))

  phecode_definitions <- phecode_definition %>%
    dplyr::filter(.data$phecode %in% unique(phecodes_mapped$phecode)) %>%
    dplyr::filter(.data$phecode %in% unique(manifest_phecodes$phecode))

  # Phecode phenotype function -----------------------------------------------------
  # these are the variables that are used in the mapply function
  phewas_ID <- phecode_definitions[["phecode_name"]]
  code_of_interest <- phecode_definitions[["phecode"]]
  l_1 <- phecode_definitions[["l1"]]
  u_1 <- phecode_definitions[["u1"]]
  l_2 <- phecode_definitions[["l2"]]
  u_2 <- phecode_definitions[["u2"]]
  l_3 <- phecode_definitions[["l3"]]
  u_3 <- phecode_definitions[["u3"]]
  sex <- phecode_definitions[["sex_num"]]

  if (isFALSE(no_phecodes)) {
    if (is.numeric(N_cores)) {
      epi_phenotypes <- parallel::mcmapply(phecoding,
                                           phewas_ID,code_of_interest,l_1,u_1,l_2,u_2,l_3,u_3,sex,
                                           MoreArgs = list(j=phecodes_mapped, k=control_exclusions, l=combined_sex),
                                           SIMPLIFY = FALSE,USE.NAMES = T,mc.cores = N_cores)
    } else {
      epi_phenotypes <- mapply(phecoding,
                                     phewas_ID,code_of_interest,l_1,u_1,l_2,u_2,l_3,u_3,sex,
                                     MoreArgs = list(j=phecodes_mapped, k=control_exclusions, l=combined_sex),
                                     SIMPLIFY = FALSE,USE.NAMES = T)
    }
    # save
    saveRDS(epi_phenotypes,phecode_save_location)
  }
  # Exclusion range function ------------------------------------------------
  # Now making a file that contains just the 233 groups of exclusions. This is then used in the primary care derived curated phenotypes where appropriate
  if (isFALSE(no_range_ID)) {
    # this uses the previously defined phecode definitions as input
    exclude_definitions <- phecode_definitions %>%
      dplyr::group_by(.data$range_ID) %>%
      dplyr::summarise(l1=max(.data$l1),u1=max(.data$u1),l2=max(.data$l2),u2=max(.data$u2),l3=max(.data$l3),u3=max(.data$u3))
    l1_e <- exclude_definitions[["l1"]]
    u1_e <- exclude_definitions[["u1"]]
    l2_e <- exclude_definitions[["l2"]]
    u2_e <- exclude_definitions[["u2"]]
    l3_e <- exclude_definitions[["l3"]]
    u3_e <- exclude_definitions[["u3"]]
    names_e <-exclude_definitions[["range_ID"]]
    # run
    range_exclusions <- mapply(exclusion_generator,
                                     names_e,
                                     l1_e,
                                     u1_e,
                                     l2_e,
                                     u2_e,
                                     l3_e,
                                     u3_e,
                                     MoreArgs = list(h=phecodes_mapped),
                                     SIMPLIFY = F, USE.NAMES = T)
    # create an empty exclusion range
    eid <- ""
    R_0 <- data.frame(eid)
    # append for final output saved as list object
    range_exclusions_edit <- append(range_exclusions,list(R_0=R_0))
    saveRDS(range_exclusions_edit,range_ID_save_location)
  }
}
