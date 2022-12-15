#' Makes primary care quantitative phenotypes guided by PheWAS manifest
#'
#' @param a PheWAS_ID
#' @param b primary care code list
#' @param c existence of limits or not
#' @param d lower_limit quant value
#' @param e upper_limit quant value
#' @param prim_care primary care data
#' @param DOB DOB information
#' @param code_lists Full file path of the folder containing alternative primary care code lists.
#' @return A list of primary care quantitative phenotypes
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data

quant_primary_care <- function(a,b,c,d,e,prim_care,DOB,code_lists) {
# message to track progress
message(a)
# create an age column for that phenotype
PheWAS_ID_age <- paste0(a,"_age")
# read in codes from list
if(code_lists=="default"){
codes <-  data.table::fread(system.file("extdata","PQP_codes",b, package = "DeepPheWAS"))
} else {
  codes <-  data.table::fread(paste0(code_lists,"/",b))
}

if(c==0) {
  # splitting by Read code Version
  read_V2 <- codes %>%
    dplyr::filter(.data$type=="read_V2") %>%
    dplyr::pull(.data$read_code)
  read_V3 <- codes %>%
    dplyr::filter(.data$type=="read_V3") %>%
    dplyr::pull(.data$read_code)
  # extracting values
  prim_care_values <- prim_care %>%
    dplyr::filter(.data$read_2 %in% read_V2 | .data$read_3 %in% read_V3) %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::filter(.data$value1 != "",
                  !grepl("^OPR|\\^",.data$value1),
                  !grepl("^OPR|\\^",.data$value2),
                  !grepl("^OPR|\\^",.data$value3)) %>%
    dplyr::mutate(date=.data$event_dt,
                  value1=as.numeric(.data$value1)) %>%
    dplyr::arrange(.data$date) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(DOB) %>%
    dplyr::mutate(age=round(lubridate::time_length(x = difftime(.data$date,.data$DOB),unit = "years"),digits = 0)) %>%
    dplyr::select(.data$eid,.data$value1,.data$date,.data$age) %>%
    purrr::set_names(c("eid",a,"earliest_date",PheWAS_ID_age))
  return(prim_care_values)
} else if (c==1) {
  # splitting by Read code Version
  read_V2 <- codes %>%
    dplyr::filter(.data$type=="read_V2") %>%
    dplyr::pull(.data$read_code)
  read_V3 <- codes %>%
    dplyr::filter(.data$type=="read_V3") %>%
    dplyr::pull(.data$read_code)
  # extracting values
  prim_care_values <- prim_care %>%
    dplyr::filter(.data$read_2 %in% read_V2 | .data$read_3 %in% read_V3) %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::filter(.data$value1 != "",
                  !grepl("^OPR|\\^",.data$value1),
                  !grepl("^OPR|\\^",.data$value2),
                  !grepl("^OPR|\\^",.data$value3)) %>%
    dplyr::mutate(date=.data$event_dt,
                  value1=as.numeric(.data$value1)) %>%
    dplyr::filter(.data$value1>d,.data$value1<e) %>%
    dplyr::arrange(.data$date) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(DOB) %>%
    dplyr::mutate(age=round(lubridate::time_length(x = difftime(.data$date,.data$DOB),unit = "years"),digits = 0)) %>%
    dplyr::select(.data$eid,.data$value1,.data$date,.data$age) %>%
    purrr::set_names(c("eid",a,"earliest_date",PheWAS_ID_age))
  return(prim_care_values)
}
}

#' MMakes primary care quantitative phenotypes guided by PheWAS manifest.
#'
#' @param GPC Full path of the primary care clinical data file from UK Biobank.
#' @param DOB Full path of a file describing participant date of birth. This does not need to be the exact date. By default for UK Biobank data, the script creates the date of birth from month and year of birth (data-fields 52 and 34).
#' @param phenotype_save_file Full path of the save file for the generated concepts RDS to be used for phenotype creation.
#' @param N_cores Number of cores requested if parallel computing is desired. Defaults to single core computing.
#' @param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#' @param code_list_overide Full file path of the folder containing alternative primary care code lists.
#' @return An RDS file containing lists of data-frames each representing a single primary-care quantitative phenotype
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

primarycare_quantitative_phenotypes <- function(GPC,
                                                DOB,
                                                phenotype_save_file,
                                                N_cores,
                                                PheWAS_manifest_overide,
                                                code_list_overide){
# Load in defining variables
# primary care data
  if(!base::file.exists(GPC)){
    rlang::abort(base::paste0("'GPC' must be a file"))
  }

  GP_C <- data.table::fread(GPC)
  GPC_expected_colname <- c("eid","data_provider","event_dt","read_2","read_3","value1","value2","value3")

  if(!all(tibble::has_name(GP_C,GPC_expected_colname))){
    base::warning(base::paste0("'GPC' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),base::paste(GPC_expected_colname, collapse=","),paste0("'
                                 not:
                                 "),base::paste(colnames(GP_C), collapse=","),
                  base::paste0("
                                 differences between inputed file and expected are:
                                 "),base::paste(setdiff_all(names(GP_C),GPC_expected_colname), collapse=","))

  }
  prim_care <- GP_C

# DOB
if(!file.exists(DOB)){
  rlang::abort(paste0("'DOB' must be a file"))
}
DOB_df <- data.table::fread(DOB)
DOB_df_colname <- c("eid","DOB")

if(!all(tibble::has_name(DOB_df,DOB_df_colname))){
  warning(paste0("'DOB' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(DOB_df_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(DOB_df), collapse=","),
          paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(DOB_df),DOB_df_colname), collapse=","))

}

DOB <- DOB_df

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
# N_Cores
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
# save location
phenotype_save_location <- phenotype_save_file
new_folder <- stringr::str_remove(phenotype_save_file,basename(phenotype_save_file))
if(!dir.exists(new_folder)){
  dir.create(new_folder)}

# selects phenotypes based on the phenotype manifest.
primary_care_phewas_ID_info <- PheWAS_manifest %>%
  dplyr::filter(.data$category=="primary_care_quantitative_phenotype") %>%
  dplyr::mutate(primary_care_code_list=paste0(.data$primary_care_code_list,".gz"))
# Run functions
if(is.numeric(N_cores)) {
  primary_care_quantitiative_data <- parallel::mcmapply(quant_primary_care,
                                                        primary_care_phewas_ID_info$PheWAS_ID,
                                                        primary_care_phewas_ID_info$primary_care_code_list,
                                                        primary_care_phewas_ID_info$limits,
                                                        primary_care_phewas_ID_info$lower_limit,
                                                        primary_care_phewas_ID_info$upper_limit,
                                                        MoreArgs = list(prim_care=prim_care,DOB=DOB,code_lists=code_list_overide),
                                                        SIMPLIFY = F,mc.cores = N_cores,USE.NAMES = T)
} else {
  primary_care_quantitiative_data <- mapply(quant_primary_care,
                                            primary_care_phewas_ID_info$PheWAS_ID,
                                            primary_care_phewas_ID_info$primary_care_code_list,
                                            primary_care_phewas_ID_info$limits,
                                            primary_care_phewas_ID_info$lower_limit,
                                            primary_care_phewas_ID_info$upper_limit,
                                            MoreArgs = list(prim_care=prim_care,DOB=DOB,code_lists=code_list_overide),
                                            SIMPLIFY = F,USE.NAMES = T)

}
saveRDS(primary_care_quantitiative_data,phenotype_save_location)
}
