#' Extracts codes for concepts into the various versions required for other functions
#' @param a source of data
#' @param code_list list of codes
#' @param health_data dataframe containing health_data
#' @return A dataframe containing of extracted data, each data item is a row with participant ID, code and source
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=

concept_source_data_extraction <- function(a,code_list,health_data) {
  . <- NULL

  all_codes <- code_list %>%
    dplyr::filter(.data$Source==a) %>%
    dplyr::mutate(Code=as.character(.data$Code))

  included_codes  <- all_codes %>%
    dplyr::filter(.data$Source==a,.data$Decision==1)

  heath_data_source <- health_data %>%
    dplyr::filter(.data$source==a)

  # and looks those up in the database
  all_codes_extracted <- heath_data_source %>%
    dplyr::filter(.data$code %in% all_codes$Code)

  included_codes_extracted <- heath_data_source %>%
    dplyr::filter(.data$code %in% included_codes$Code)

  ### wide date
  min_date_label <- paste0("first_code_",a)

  wide_date <- included_codes_extracted %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::summarise(!!min_date_label:=min(.data$date))

  #### wide file
  N_source_type <- included_codes_extracted %>%
    dplyr::count(.data$eid) %>%
    dplyr::rename(!!a:=.data$n)

  #### all dates
  all_dates <- included_codes_extracted %>%
    dplyr::select(.data$eid,.data$date)

  #### TN
  included_ID <- unique(included_codes_extracted$eid)
  included_ID_N <-length(included_ID)
  type_of_code <- a
  half_TN <- data.frame(included_ID_N,type_of_code) %>%
    dplyr::select(Source=.data$type_of_code,Total_participants=.data$included_ID_N)

  #### CTID
  N_codes <- all_codes_extracted %>%
    dplyr::group_by(.data$code) %>%
    dplyr::count() %>%
    dplyr::mutate(Source={{a}},) %>%
    dplyr::rename(N_Codes=.data$n)

  N_codes_ID <- all_codes_extracted %>%
    dplyr::group_by(.data$code) %>%
    dplyr::summarise(N_Participants= dplyr::n_distinct(.data$eid))

  # and identifies how many participants have only one code out of the code list in their records
  IDs_only_one_code <- all_codes_extracted[all_codes_extracted$eid %in% names(which(table(all_codes_extracted$eid) < 2)), ]
  only_one_code <- IDs_only_one_code %>%
    dplyr::group_by(.data$code) %>%
    dplyr::summarise(N_Participants_with_1_code= dplyr::n_distinct(.data$eid))

  CTID_primer <- all_codes %>%
    dplyr::left_join(N_codes,by=c("Code"="code","Source")) %>%
    dplyr::left_join(N_codes_ID,by=c("Code"="code")) %>%
    dplyr::left_join(only_one_code,by=c("Code"="code")) %>%
    dplyr::select(.data$Code,.data$Description,.data$Source,.data$Decision,.data$N_Codes,.data$N_Participants,.data$N_Participants_with_1_code)

  output <- list(wide_date=wide_date,wide_file=N_source_type,all_dates=all_dates,half_TN=half_TN,included_IDs=included_ID,CTID_primer=CTID_primer)

  return(output)

}

#' Minor function for concept creation
#'
#' @param x numeric list
#' @param included_ID list of codes
#' @return A dataframe containing of extracted data, each data item is a row with participant ID, code and source
#' @keywords internal

calc_minus_source_ID <- function(x,included_ID) {
  current_ID <- unlist(included_ID[x])
  alternate_ID <- unlist(included_ID[-x])
  exclusive_ID <- sum((!current_ID %in% alternate_ID),na.rm = T)

  return(exclusive_ID)
}

#' Creates clinical concepts for combination in composite phenotypes
#'
#' @param x PheWAS_ID for concept
#' @param z location of the file
#' @param health_data dataframe containing health_data
#' @param code_list_folder marker to detect if default concept data is used or alternative folder of concept code lists is provided
#' @return A list each containing 2 lists one for WDA and one for all dates used later in the pipeline for clinical concepts
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=
clinical_code_lookup <- function (x,z,health_data,code_list_folder) {
  . <- NULL
  message(x)
  # this is the rated code list for each concept
  if(code_list_folder=="default"){
    code_list_alpha <- data.table::fread(system.file("extdata","concept_codes",z, package = "DeepPheWAS"))
  } else {
    code_list_alpha <- data.table::fread(z)
  }
  concept_name <- x
  sources_alpha <- unique(code_list_alpha$Source)
  # variable if statement to isolate presence or absence of ICD 10/9 or both
  if(any(stringr::str_detect("ICD10$",sources_alpha),na.rm = T) & any(stringr::str_detect("ICD9$",sources_alpha),na.rm = T)){
    # filter to remove ICD10 or 9 or both
    no_ICD <- code_list_alpha %>%
      dplyr::filter(.data$Source!="ICD10", .data$Source!="ICD9")
    # create a seperate data frame for each position of ICD code in 10/9 or both
    ICD10_1 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_1")
    ICD10_2 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_2")
    ICD10_3 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_3")
    ICD9_1 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_1")
    ICD9_2 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_2")
    ICD9_3 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_3")
    # combine into a new code_list variable that is used throughout function
    code_list <- no_ICD %>%
      dplyr::bind_rows(ICD10_1,ICD10_2,ICD10_3,ICD9_1,ICD9_2,ICD9_3)
    sources <- unique(code_list$Source)
  } else if(any(stringr::str_detect("ICD10$",sources_alpha),na.rm = T)) {
    no_ICD <- code_list_alpha %>%
      dplyr::filter(.data$Source!="ICD10")
    ICD10_1 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_1")
    ICD10_2 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_2")
    ICD10_3 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD10") %>%
      dplyr::mutate(Source="ICD10_3")
    code_list <- no_ICD %>%
      dplyr::bind_rows(ICD10_1,ICD10_2,ICD10_3)
    sources <- unique(code_list$Source)
  } else if(any(stringr::str_detect("ICD9$",sources_alpha),na.rm = T)){
    no_ICD <- code_list_alpha %>%
      dplyr::filter(.data$Source!="ICD9")
    ICD9_1 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_1")
    ICD9_2 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_2")
    ICD9_3 <- code_list_alpha %>%
      dplyr::filter(.data$Source=="ICD9") %>%
      dplyr::mutate(Source="ICD9_3")
    code_list <- no_ICD %>%
      dplyr::bind_rows(ICD9_1,ICD9_2,ICD9_3)
    sources <- unique(code_list$Source)
  } else {
    code_list <- code_list_alpha
    sources <-  sources_alpha
  }

  per_source_data <- mapply(concept_source_data_extraction,
                            sources,
                            MoreArgs = list(code_list=code_list,health_data=health_data),
                            SIMPLIFY = F)

  rowMins <- function(x, na.rm = FALSE) {
    stopifnot(is.data.frame(x))
    do.call(pmin, c(unname(as.list(x)), list(na.rm = na.rm)))
  }

  ## Wide dates all file primary file used for phenotype generation
  wide_dates <- lapply(per_source_data,'[[',"wide_date") %>%
    purrr::reduce(dplyr::full_join) %>%
    dplyr::mutate(earliest_date = rowMins(dplyr::across(2:dplyr::last_col()), na.rm = TRUE)) %>%
    dplyr::select(.data$eid,.data$earliest_date)

  wide_count <- lapply(per_source_data,'[[',"wide_file") %>%
    purrr::reduce(dplyr::full_join) %>%
    dplyr::mutate(any_code = rowSums(dplyr::across(2:dplyr::last_col()), na.rm = TRUE))

  WA <- wide_dates %>%
    dplyr::left_join(wide_count)

  WA[is.na(WA)] <- 0

  WDA <- WA %>%
    dplyr::select(.data$eid,!!x:=.data$any_code,.data$earliest_date)

  ##  all dates
  AD <- lapply(per_source_data,'[[',"all_dates") %>%
    purrr::reduce(dplyr::bind_rows)

  ## AID
  AID <- data.frame(length(unique(lapply(per_source_data,'[[',"included_IDs") %>%
                                    purrr::reduce(c)))) %>%
    dplyr::rename(N_participants=1)

  ##  TN
  half_TN <- lapply(per_source_data,'[[',"half_TN") %>%
    purrr::reduce(dplyr::bind_rows)

  included_ID <- lapply(per_source_data,'[[',"included_IDs")

  range_of_x <- c(1:length(included_ID))
  only_source_ID <- mapply(calc_minus_source_ID,range_of_x,MoreArgs = list(included_ID=included_ID),SIMPLIFY = F) %>%
    purrr::reduce(c) %>%
    as.data.frame() %>%
    dplyr::rename(N_Participants_only_this_code_source=1)

  TN <- half_TN %>%
    dplyr::bind_cols(only_source_ID)

  ## CTID
  CTID_primer <- lapply(per_source_data,'[[',"CTID_primer") %>%
    purrr::reduce(dplyr::bind_rows)
  CTID_primer[is.na(CTID_primer)] <- 0

  CTID <- CTID_primer %>%
    dplyr::arrange(dplyr::desc(.data$Decision),dplyr::desc(.data$N_Participants),.data$Source)

  #  fwrite(WA,paste0(concept_save_location,"/",paste0(concept_name,"_WA.txt")), sep = "\t")
  #  fwrite(WDA,paste0(concept_save_location,"/",paste0(concept_name,"_WDA.txt")),quote = T, na = NA, sep = "\t")
  #  fwrite(AID,paste0(concept_save_location,"/",paste0(concept_name,"_AID.txt")), sep = "\t")
  #  fwrite(AD,paste0(concept_save_location,"/",paste0(concept_name,"_AD.txt")), sep = "\t")
  #  fwrite(CTID,paste0(concept_save_location,"/",paste0(concept_name,"_CTID.txt")),quote = T, na = NA, sep = "\t")

  new_return <- list(WDA=WDA,AD=AD)
  return(new_return)
}

#' Creates prescription concepts for combination in composite phenotypes
#'
#' @param x PheWAS_ID for concept
#' @param y Read_V2 drug list
#' @param z Drug key search terms
#' @param GP_P dataframe containing GP prescription data
#' @param V2_drugs dataframe containing information on drugs coded using Read_V2 codes.
#' @param code_list_folder marker to detect if default concept data is used or alternative folder of concept code lists is provided
#' @return A list each containing 2 lists one for WDA and one for all dates used later in the pipeline for prescription concepts
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=
drug_code_lookup <- function(x,y,z,GP_P,V2_drugs,code_list_folder) {
  . <- NULL
  message(x)
  ## creating a name for the files
  both_names <- x

  if(code_list_folder=="default"){
    drug_list_BNF_DMD <- data.table::fread(system.file("extdata","concept_codes",z, package = "DeepPheWAS"),sep = ",") %>%
      dplyr::mutate(drug_unique=tolower(.data$drug_unique))%>%
      dplyr::pull()

    drug_df_BNF_DMD <- data.table::fread(system.file("extdata","concept_codes",z, package = "DeepPheWAS"),sep = ",") %>%
      dplyr::mutate(drug_unique=tolower(.data$drug_unique))

    drug_list_V2 <- data.table::fread(system.file("extdata","concept_codes",y, package = "DeepPheWAS")) %>%
      dplyr::filter(.data$Decision==1) %>%
      dplyr::select(1) %>%
      dplyr::pull()
  } else {
    drug_list_BNF_DMD <- data.table::fread(z,sep = ",") %>%
      dplyr::mutate(drug_unique=tolower(.data$drug_unique)) %>%
      dplyr::pull()

    drug_df_BNF_DMD <- data.table::fread(z,sep = ",") %>%
      dplyr::mutate(drug_unique=tolower(.data$drug_unique))

    drug_list_V2 <- data.table::fread(y) %>%
      dplyr::filter(.data$Decision==1) %>%
      dplyr::select(1) %>%
      dplyr::pull()
  }

 ## use these names in a string search, this not only searches for the keywords but also reports which keywords are then found in each row. I use this
  ## to add a keyword for describing the number of participants found for each broad keyword.
    drug_text_search <- GP_P %>%
    dplyr::mutate(drug_name=tolower(.data$drug_name)) %>%
    dplyr::filter(grepl(paste(drug_list_BNF_DMD, collapse="|"),.data$drug_name)) %>%
    dplyr::mutate(keyword_or_V2code=stringr::str_extract(.data$drug_name,paste(drug_list_BNF_DMD, collapse="|")),
                  keyword_or_V2code=as.character(.data$keyword_or_V2code),
                  keyword_or_V2code=stringr::str_to_title(.data$keyword_or_V2code),
                  Source="BNF/DMD")

  ## now use read_v2 codes to find the small percentage of codes with only read_v2 identification this only searches in those with read_v2 codes and not those eids
  ## already found within the keyword search above
  drug_V2_search <- GP_P %>%
    dplyr::filter(.data$read_2!="NA", !.data$row_number %in% drug_text_search$row_number, .data$read_2 %in% drug_list_V2) %>%
    dplyr::mutate(keyword_or_V2code=.data$read_2,
                  Source="Read_V2") %>%
    dplyr::left_join(V2_drugs, by=c("read_2"="read_code")) %>%
    dplyr::select(-.data$status_flag)

  ## combine for a list of results
  all_drug_results <- dplyr::bind_rows(drug_text_search,drug_V2_search) %>%
    dplyr::distinct(.data$row_number, .keep_all = TRUE)

  ## create all_id
  all_id_drugs <- list(length(unique(all_drug_results$eid)))

  ## This creates the codes total ID's or CTID files for prescription data, it does this for the keyword/Read_v2 code
  ## first it counts all incidents of the codes
  CTID_ish <- all_drug_results %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::count() %>%
    dplyr::rename(N_Codes=.data$n)

  ## then counts the N of unique eids i.e. how many particpants in total this represents
  CTID2_ish <- all_drug_results %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::summarise(N_Participants=dplyr::n_distinct(.data$eid))

  ## identifies those with only one code in the record
  prescription_only_one_code <- all_drug_results %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::summarise(N_eid=dplyr::n()) %>%
    dplyr::filter(.data$N_eid <2)

  ## and uses the info above to show how many for each keyword/Read_v2 code have only a single code
  prescription_ids_ooc <- all_drug_results %>%
    dplyr::filter(.data$eid %in% prescription_only_one_code$eid) %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::summarise(N_Particicpants_with_1_code= dplyr::n_distinct(.data$eid))

  ## calculates how many prescriptions for each code there is on average
  CTID3_ish <- all_drug_results %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::summarise(Mean_prescriptions_per_participant=round((dplyr::n()/dplyr::n_distinct(.data$eid)),digits = 1))

  ## this does the same as the above but removing those with only one prescription in their records
  CTID4_ish <- all_drug_results %>%
    dplyr::filter(!.data$eid %in% prescription_only_one_code$eid) %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::summarise(`Mean_prescriptions_per_participant*`=round((dplyr::n()/dplyr::n_distinct(.data$eid)),digits=1))

  ## the two below add in data on description and Source for the final table
  read_descriptions <- drug_V2_search %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::summarise(terms_x=unique(.data$term_description),Source=unique(.data$Source))

  drug_list_description <- drug_df_BNF_DMD %>%
    dplyr::mutate(Source="BNF/DMD",
                  Description=.data$drug_unique) %>%
    dplyr::rename(keyword_or_V2code=1)

  ## making the final CTID table
  CTID_drugs <- CTID_ish %>%
    dplyr::left_join(CTID2_ish)  %>%
    dplyr::left_join(prescription_ids_ooc) %>%
    dplyr::left_join(CTID3_ish) %>%
    dplyr::left_join(CTID4_ish)  %>%
    dplyr::left_join(read_descriptions) %>%
    dplyr::left_join(drug_list_description, by="keyword_or_V2code") %>%
    dplyr::mutate(Source=ifelse(is.na(.data$Source.x),.data$Source.y,.data$Source.x),
                  Description=ifelse(is.na(.data$Description),.data$terms_x,.data$Description),
                  Decision=1) %>%
    dplyr::select(1,.data$Description,.data$Source,.data$Decision,.data$N_Codes,.data$N_Participants,.data$N_Particicpants_with_1_code,.data$Mean_prescriptions_per_participant,.data$`Mean_prescriptions_per_participant*`) %>%
    dplyr::mutate_all(~replace(., is.na(.), 0))

  ## now formatting dates
  all_drugs_with_date <- all_drug_results %>%
    dplyr::filter(.data$issue_date!="NA") %>%
    dplyr::rename(code_date=.data$issue_date) %>%
    dplyr::mutate(code_date=lubridate::ymd(.data$code_date))
  all_drugs_with_date <- all_drugs_with_date %>%
    tidyr::drop_na(.data$code_date)

  ## create earliest and latest date of prescription
  all_drugs_time <- all_drugs_with_date %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::summarise(earliest_date=min(.data$code_date), any_code=max(.data$code_date), any_code=dplyr::n())

  ## all dates file
  all_dates_drugs <- all_drugs_with_date %>%
    dplyr::select(.data$eid,.data$code_date)

  ## to make equivalent of all_wide_dates and all_wide file need to create a separation in the codes into BNF/DMD and Read_V2
  ## do for BNF
  BNF_DMD <- all_drug_results %>%
    dplyr::filter(.data$Source=="BNF/DMD")
  BNF_DMD_wide <- BNF_DMD %>%
    dplyr::count(.data$eid,.data$keyword_or_V2code) %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_wider(names_from = .data$keyword_or_V2code, values_from = .data$n) %>%
    dplyr::select(-.data$row) %>%
    dplyr::group_by(.data$eid)

  BNF_DMD_wide[is.na(BNF_DMD_wide)] <- 0

  BNF_DMD_wide <- BNF_DMD_wide %>%
    dplyr::summarise_all(list(max)) %>%
    dplyr::mutate(BNF_DMD = rowSums(.[-1],na.rm = TRUE))

  BNF_DMD_wide_date <- BNF_DMD %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::summarise(BNF_DMD_first_code=min(.data$issue_date))

  ## Now for Read_V2
  Read_V2 <- all_drug_results %>%
    dplyr::filter(.data$Source=="Read_V2")

  Read_V2_wide <- Read_V2 %>%
    dplyr::count(.data$eid,.data$keyword_or_V2code) %>%
    dplyr::group_by(.data$keyword_or_V2code) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_wider(names_from = .data$keyword_or_V2code, values_from = .data$n) %>%
    dplyr::select(-.data$row) %>%
    dplyr::group_by(.data$eid)

  Read_V2_wide[is.na(Read_V2_wide)] <- 0

  Read_V2_wide <- Read_V2_wide %>%
    dplyr::summarise_all(list(max)) %>%
    dplyr::mutate(Read_V2 = rowSums(.[-1],na.rm = TRUE))

  Read_V2_wide_date <- Read_V2 %>%
    dplyr::group_by(.data$eid) %>%
    dplyr::summarise(Read_V2_first_code=min(.data$issue_date))

  ## Now combine them
  wide_drugs <- BNF_DMD_wide %>%
    dplyr::left_join(Read_V2_wide) %>%
    dplyr::select(.data$eid,.data$BNF_DMD,.data$Read_V2)

  wide_drugs[is.na(wide_drugs)] <- 0

  wide_drugs <- wide_drugs %>%
    dplyr::mutate(any_code = .data$BNF_DMD + .data$Read_V2)

  ##now with dates added in
  wide_drugs_date <- wide_drugs %>%
    dplyr::left_join(BNF_DMD_wide_date) %>%
    dplyr::left_join(Read_V2_wide_date) %>%
    dplyr::mutate(earliest_date=ifelse(.data$BNF_DMD_first_code==0, .data$Read_V2_first_code, .data$BNF_DMD_first_code)) %>%
    dplyr::select(.data$eid,.data$earliest_date,.data$BNF_DMD,.data$Read_V2,.data$any_code) %>%
    dplyr::filter(.data$earliest_date!="NA")

  wide_dates_all <- wide_drugs_date %>%
    dplyr::select(.data$eid,!!x:=.data$any_code,.data$earliest_date)

  #Total Numbers
  BNF_DMD_ID <- length(unique(BNF_DMD$eid))
  Read_V2_ID <- length(unique(Read_V2$eid))

  ID_N <- c(BNF_DMD_ID,Read_V2_ID)
  just_ID_N <- c(BNF_DMD_ID,Read_V2_ID)
  numbers_row <- c("BNF/DMD","Read_V2")

  total_numbers_drugs <- data.table::setnames(data.frame(matrix(ncol = 1, nrow = 2)), "flib") %>%
    tibble::add_column(ID_N,just_ID_N,numbers_row) %>%
    dplyr::select(Source=.data$numbers_row,Total_Participants=.data$ID_N,N_Participants_only_this_code_source=.data$just_ID_N)

  #  fwrite(wide_drugs_date,paste0(concept_save_location,"/",paste0(both_names,"_WA.txt")), quote = TRUE, na=NA)
  #  fwrite(all_dates_drugs,paste0(concept_save_location,"/",paste0(both_names,"_AD.txt")), quote = TRUE, na=NA)
  #  fwrite(CTID_drugs,paste0(concept_save_location,"/",paste0(both_names,"_CTID.txt")), quote = TRUE, na=NA)
  #  fwrite(all_id_drugs,paste0(concept_save_location,"/",paste0(both_names,"_AID.txt")), quote = TRUE, na=NA)
  #  fwrite(total_numbers_drugs,paste0(concept_save_location,"/",paste0(both_names,"_TN.txt")), quote = TRUE, na=NA)
  #  fwrite(wide_dates_all,paste0(concept_save_location,"/",paste0(both_names,"_WDA.txt")), quote = TRUE, na=NA)

  new_return <- list(WDA=wide_dates_all,AD=all_dates_drugs)

  return(new_return)
}

#' Creates prescription concepts for combination in composite phenotypes
#'
#' @param health_data Full file path of the min_data file. Selecting will generate clinical concepts.
#' @param GPP Full file path of the GP prescription data. Selecting will generate prescription concepts.
#' @param concept_save_file Full file path for the save file for the generated concepts RDS used for concept creation.
#' @param all_dates_save_file Full file path for the save file for the generated all_dates.RDS used for per-event combinations of concepts.
#' @param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#' @param code_list_folder_override Full file for the folder containing code lists only use if not using default stored in R package.
#' @return Two RDS objects for concepts one containing all_dates and one containing summarised N codes and earliest date.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

creating_concepts <- function(health_data,
                              GPP,
                              concept_save_file,
                              all_dates_save_file,
                              PheWAS_manifest_overide,
                              code_list_folder_override){
concepts <- list()
all_dates <- list()
# defining save locations for
# concept.RDS

phenotype_save_location <- concept_save_file
new_folder <- stringr::str_remove(concept_save_file,basename(concept_save_file))
if(!dir.exists(new_folder)){
  dir.create(new_folder)}

# all dates
all_dates_save_location <- all_dates_save_file
new_folder <- stringr::str_remove(all_dates_save_file,basename(all_dates_save_file))
if(!dir.exists(new_folder)){
  dir.create(new_folder)}

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

# code list folder location
if(is.null(code_list_folder_override)){
  code_list_folder <- "default"
} else {
  code_list_folder <- "new"
  code_list_folder_location <- code_list_folder_override
  if(!dir.exists(code_list_folder_location)){
    rlang::abort(paste0("'code_list_folder_override' must be an active directory"))
  }
}

if(!is.null(health_data)) {
  # health_data
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

  # Clinical search codes works by searching for all clinical code list files in specified folder
  if(code_list_folder=="default"){
    concept_file_list <- data.table::fread(system.file("extdata","clinical_concept_file_list.csv.gz", package = "DeepPheWAS")) %>%
      dplyr::pull()
  } else {
    concept_file_list <- list.files(code_list_folder_location,pattern = "_codes_rated.csv.gz")
  }

  concept_clinical_codes <- data.frame(names=(stringr::str_remove(concept_file_list,"_codes_rated.csv.gz")))%>%
    dplyr::left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%
    tidyr::drop_na(.data$PheWAS_ID) %>%
    dplyr::select(.data$names) %>%
    dplyr::pull()

  # gets PheWAS_ID which is used to save and then access the files later from joining with the PheWAS_manifest
  concept_clinical_codes_PheWAS_ID <- data.frame(names=(stringr::str_remove(concept_file_list,"_codes_rated.csv.gz"))) %>%
    dplyr::left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%
    tidyr::drop_na(.data$PheWAS_ID) %>%
    dplyr::select(.data$PheWAS_ID) %>%
    dplyr::pull()

  # running the clinical concept function
  if(code_list_folder=="default"){
    concept_clinical_codes_list <- lapply(concept_clinical_codes, function(x) paste0(x,"_codes_rated.csv.gz"))
  } else {
    concept_clinical_codes_list <- lapply(concept_clinical_codes, function(x) paste0(code_list_folder_location,"/",x,"_codes_rated.csv.gz"))
  }
  clinical_concepts <- mapply(clinical_code_lookup,
                              concept_clinical_codes_PheWAS_ID,
                              concept_clinical_codes_list,
                              MoreArgs = list(health_data=health_data,code_list_folder=code_list_folder),
                              SIMPLIFY = F,
                              USE.NAMES = T)
  # is a list with two outputs, retrive WDA which saves as concepts.RDS and AD which saves as all_dates.RDS
  WDA <- lapply(clinical_concepts,'[[',"WDA")
  AD <- lapply(clinical_concepts,'[[',"AD")
  new_names <- paste0(names(AD),"_AD")
  names(AD) <- new_names
  # appending lists
  concepts <- append(concepts,WDA)
  all_dates <-append(all_dates,AD)
  # for running individual code list lookups
  #l <- concept_clinical_codes[5]
  #m <- concept_clinical_codes_list[5]
  #mapply(clinical_code_lookup,l,m)
}

if(!is.null(GPP)) {
  # prescription data
  if(!file.exists(GPP)){
    rlang::abort(paste0("'GPP' must be a file"))
  }
  GPP <- data.table::fread(GPP, na.strings = "")
  GPP_expected_colname <- c("eid","data_provider","issue_date","read_2","bnf_code","dmd_code","drug_name","quantity")

  if(!all(tibble::has_name(GPP,GPP_expected_colname))){
    warning(paste0("'GPP' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(GPP_expected_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(GPP), collapse=","),
                  paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(GPP),GPP_expected_colname), collapse=","))

  }

  GP_P <- GPP %>%
    dplyr::select(.data$eid,.data$read_2,.data$issue_date,.data$drug_name) %>%
    dplyr::mutate(row_number = 1:dplyr::n())
  V2_drugs <- data.table::fread(system.file("extdata","V2_drugs.csv.gz", package = "DeepPheWAS"))
  # Prescription search codes works by searching for all clinical code list files in specified folder
  if(code_list_folder=="default"){
    prescription_concept_file_list <- data.table::fread(system.file("extdata","prescription_concept_file_list.csv.gz", package = "DeepPheWAS")) %>%
      dplyr::pull()
  } else {
    prescription_concept_file_list <- list.files(code_list_folder_location,pattern = "_BNF_DMD.csv.gz")
  }

  prescription_search_terms <- data.frame(names=(stringr::str_remove(prescription_concept_file_list,"_BNF_DMD.csv.gz"))) %>%
    dplyr::left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%
    tidyr::drop_na(.data$PheWAS_ID) %>%
    dplyr::select(.data$names) %>%
    dplyr::pull()
  # and extracting the PheWAS ID by joining with the PheWAS_manifest
  prescription_search_terms_PheWAS_ID <- data.frame(names=(stringr::str_remove(prescription_concept_file_list,"_BNF_DMD.csv.gz"))) %>%
    dplyr::left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%
    tidyr::drop_na(.data$PheWAS_ID) %>%
    dplyr::select(.data$PheWAS_ID) %>%
    dplyr::pull()

  # V2 files here unique quirk of the UK Biobank data
  if(code_list_folder=="default"){
    V2_concept_file_list <- data.table::fread(system.file("extdata","V2_concept_file_list.csv.gz", package = "DeepPheWAS")) %>%
      dplyr::pull()
  } else {
    V2_concept_file_list <- list.files(code_list_folder_location,pattern = "_V2_rated.csv.gz")
  }
  V2_search_terms <- data.frame(names=(stringr::str_remove(V2_concept_file_list,"_V2_rated.csv.gz"))) %>%
    dplyr::left_join(PheWAS_manifest, by=c("names"="concept_name")) %>%
    tidyr::drop_na(.data$PheWAS_ID) %>%
    dplyr::select(.data$names) %>%
    dplyr::pull()
  if(code_list_folder=="default"){
    prescription_search_terms_list <- lapply(prescription_search_terms, function(x) paste0(x,"_BNF_DMD.csv.gz"))
    V2_search_terms_list <- lapply(V2_search_terms, function(x) paste0(x,"_V2_rated.csv.gz"))
  } else {
    prescription_search_terms_list <- lapply(prescription_search_terms, function(x) paste0(code_list_folder_location,"/",x,"_BNF_DMD.csv.gz"))
    V2_search_terms_list <- lapply(V2_search_terms, function(x) paste0(code_list_folder_location,"/",x,"_V2_rated.csv.gz"))
  }

  # run the function
  drug_concepts <- mapply(drug_code_lookup,
                          prescription_search_terms_PheWAS_ID,
                          V2_search_terms_list,
                          prescription_search_terms_list,
                          MoreArgs = list(GP_P=GP_P,V2_drugs=V2_drugs,code_list_folder=code_list_folder),
                          SIMPLIFY = F,
                          USE.NAMES = T)
  # has two outputs as list need to extract and then append to lists
  WDA <- lapply(drug_concepts,'[[',"WDA")
  AD <- lapply(drug_concepts,'[[',"AD")
  new_names <- paste0(names(AD),"_AD")
  names(AD) <- new_names
  # append to lists
  concepts <- append(concepts,WDA)
  all_dates <-append(all_dates,AD)
}
# Save RDS ----------------------------------------------------------------
saveRDS(concepts,phenotype_save_location)
saveRDS(all_dates,all_dates_save_location)
}
