#' Extracts data-field phenotypes from study data guided by PheWAS manifest
#'
#' @param a PheWAS_ID
#' @param b QC_flag_ID
#' @param c search_id
#' @param d type of analysis
#' @param e data coding for case
#' @param f data coding for controls
#' @param g existence of limits or not
#' @param h lower_limit quant value
#' @param i upper_limit quant value
#' @param j data coding for date
#' @param k data coding for age
#' @param l selecting whether first_last_max_min_value
#' @param m selecting whether to exclude_case_control_crossover
#' @param n data coding for QC vales
#' @param min_data dataframe of min_data
#' @return A dataframe containing data field phenotype
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data

data_field_extraction <- function (a,b,c,d,e,f,g,h,i,j,k,l,m,n,min_data) {
  # a slightly different method is needed when considering traits that have a separate QC flag, we need to be able to
  # remove the measurement that has been flagged and replace with an NA before we evaluate the measurements.
  . <- NULL
  message(a)
  if(l=="first"){
    position_col <- c("first_result","matching_date_position_first")
  } else if(l=="last"){
    position_col <- c("last_result","matching_date_position_last")
  } else if(l=="min"){
    position_col <- c("min_result","matching_date_position_min")
  } else if(l=="max"){
    position_col <- c("max_result","matching_date_position_max")
  }

  age_name <- paste0(a,"_age")

  if(d=="binary") {
    # extract columns for phenotype
    per_data_field <- min_data %>%
      dplyr::select(.data$eid,tidyselect::matches(c))
    # vector used later
    data_field_colnames <- colnames(per_data_field)[-1]
    # extract date column
    date_col <- min_data %>%
      dplyr::select(.data$eid,tidyselect::matches(j))
    # vector used later
    date_col_names <- colnames(date_col)[-1]
    # need to extract cases and controls separately, to achieve this need to identify unique values and set all non-cases or non-controls as NA
    all_values <- unique(as.vector(as.matrix(per_data_field[,-1])))
    case_values <- as.numeric(unlist(strsplit(e,",")))
    control_values <- as.numeric(unlist(strsplit(f,",")))
    # leave only values that represent cases
    na_values_controls <- all_values[- which(all_values %in% case_values)]
    # and controls
    na_values_cases <- all_values[- which(all_values %in% control_values)]
    # turn non-cases/controls into NA
    just_cases <- per_data_field %>%
      dplyr::mutate(dplyr::across(dplyr::everything(),~replace(., . %in% na_values_controls,NA)))
    just_controls <- per_data_field %>%
      dplyr::mutate(dplyr::across(dplyr::everything(),~replace(., . %in% na_values_cases,NA)))
    # create cases
    cases <- just_cases %>%
      tidyr::unite(col="results",2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>%
      dplyr::mutate(results=dplyr::na_if(.data$results,"")) %>%
      tidyr::drop_na(.data$results) %>%
      dplyr::mutate(results_list=purrr::transpose(dplyr::across(tidyselect::all_of(data_field_colnames))))
    # too computationally inefficient to create all variations so create only the one that is required
    if(l=="first"){
      cases_extracted <- cases %>%
        dplyr::mutate(first_result= unlist(purrr::map(.data$results_list, ~ dplyr::first(stats::na.omit(.x)))),
                      position_first=unlist(purrr::map2(.data$results_list,.data$first_result, ~match(.y,unlist(.x)))),
                      matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_first])))+1,)
    } else if(l=="last"){
      cases_extracted <- cases %>%
        dplyr::mutate(last_result= unlist(purrr::map(.data$results_list, ~ dplyr::last(stats::na.omit(.x)))),
                      position_last=unlist(purrr::map2(.data$results_list,.data$last_result, ~match(.y,unlist(.x)))),
                      matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_last])))+1,)
    } else if(l=="min"){
      cases_extracted <- cases %>%
        dplyr::mutate(min_result= unlist(purrr::map(.data$results_list, ~ min(unlist(.x),na.rm = T))),
                      position_min=unlist(purrr::map2(.data$results_list,.data$min_result, ~match(.y,unlist(.x)))),
                      matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_min])))+1,)
    } else if(l=="max"){
      cases_extracted <- cases %>%
        dplyr::mutate(max_result= unlist(purrr::map(.data$results_list, ~ max(unlist(.x),na.rm = T))),
                      position_max=unlist(purrr::map2(.data$results_list,.data$max_result, ~match(.y,unlist(.x)))),
                      matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_max])))+1,)
    }
    cases_extracted_complete <- cases_extracted %>%
      dplyr::select(.data$eid,tidyselect::all_of(position_col)) %>%
      dplyr::rename(eid=1,pheno=2,position_match=3) %>%
      dplyr::mutate(pheno=1) %>%
      tidyr::drop_na() %>%
      dplyr::left_join(date_col, by="eid") %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
      dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                    earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
      dplyr::select(.data$eid,.data$pheno,.data$earliest_date) %>%
      dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date),
                    pheno=1) %>%
      purrr::set_names(c("eid",a,"earliest_date"))
    # create controls
    controls <- just_controls %>%
      tidyr::unite(col="results",2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>%
      dplyr::mutate(results=dplyr::na_if(.data$results,"")) %>%
      tidyr::drop_na(.data$results) %>%
      dplyr::mutate(results_list=purrr::transpose(dplyr::across(tidyselect::all_of(data_field_colnames))))
    # too computationally inefficient to create all variations so create only the one that is required
    if(l=="first"){
      controls_extracted <- controls %>%
        dplyr::mutate(first_result= unlist(purrr::map(.data$results_list, ~ dplyr::first(stats::na.omit(.x)))),
                      position_first=unlist(purrr::map2(.data$results_list,.data$first_result, ~match(.y,unlist(.x)))),
                      matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_first])))+1,)
    } else if(l=="last"){
      controls_extracted <- controls %>%
        dplyr::mutate(last_result= unlist(purrr::map(.data$results_list, ~ dplyr::last(stats::na.omit(.x)))),
                      position_last=unlist(purrr::map2(.data$results_list,.data$last_result, ~match(.y,unlist(.x)))),
                      matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_last])))+1,)
    } else if(l=="min"){
      controls_extracted <- controls %>%
        dplyr::mutate(min_result= unlist(purrr::map(.data$results_list, ~ min(unlist(.x),na.rm = T))),
                      position_min=unlist(purrr::map2(.data$results_list,.data$min_result, ~match(.y,unlist(.x)))),
                      matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_min])))+1,)
    } else if(l=="max"){
      controls_extracted <- controls %>%
        dplyr::mutate(max_result= unlist(purrr::map(.data$results_list, ~ max(unlist(.x),na.rm = T))),
                      position_max=unlist(purrr::map2(.data$results_list,.data$max_result, ~match(.y,unlist(.x)))),
                      matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_max])))+1,)
    }
    controls_extracted_complete <- controls_extracted %>%
      dplyr::select(.data$eid,tidyselect::all_of(position_col)) %>%
      dplyr::rename(eid=1,pheno=2,position_match=3) %>%
      dplyr::mutate(pheno=1) %>%
      tidyr::drop_na() %>%
      dplyr::left_join(date_col, by="eid") %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
      dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                    earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
      dplyr::select(.data$eid,.data$pheno,.data$earliest_date) %>%
      dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date),
                    pheno=0) %>%
      purrr::set_names(c("eid",a,"earliest_date"))
    # union of cases and controls to ensure no repeat ID's
    union_cases_controls <- cases_extracted_complete %>%
      dplyr::inner_join(controls_extracted_complete, by="eid") %>%
      dplyr::pull(.data$eid)
    if(m==1) {
      # exclude crossover cases/controls from being either cases or controls
      phenotype_extraction_final <- cases_extracted_complete %>%
        dplyr::bind_rows(controls_extracted_complete) %>%
        dplyr::filter(!.data$eid %in% union_cases_controls)
    } else if(m==0){
      #exclude the crossover only from controls to ensure that there are no duplicate ids
      controls_modified <- controls_extracted_complete %>%
        dplyr::filter(!.data$eid %in% union_cases_controls)
      phenotype_extraction_final <- cases_extracted_complete %>%
        dplyr::bind_rows(controls_modified)
    }
    return(phenotype_extraction_final)
  } else {
    if(is.na(b) || ncol(min_data %>% dplyr::select(.data$eid,tidyselect::matches(stats::na.omit(b)))) <2){
      # extract columns for phenotype
      per_data_field <- min_data %>%
        dplyr::select(.data$eid,tidyselect::matches(c))
      # vector used later
      data_field_colnames <- colnames(per_data_field)[-1]
      # extract date column
      date_col <- min_data %>%
        dplyr::select(.data$eid,tidyselect::matches(j))
      # vector used later
      date_col_names <- colnames(date_col)[-1]
      # dplyr::selecting single value for quantitative measure where multiple options are available,
      # joins with age at assessment to get accurate age when measure taken.
      phenotype_extraction <- per_data_field %>%
        tidyr::unite(col="results",2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>%
        dplyr::mutate(results=dplyr::na_if(.data$results,"")) %>%
        tidyr::drop_na(.data$results) %>%
        dplyr::mutate(results_list = purrr::transpose(dplyr::across(tidyselect::all_of(data_field_colnames))))

      # too computationally inefficient to create all variations so create only the one that is required
      if(l=="first"){
        phenotype_extracted <- phenotype_extraction %>%
          dplyr::mutate(first_result= unlist(purrr::map(.data$results_list, ~ dplyr::first(stats::na.omit(.x)))),
                        position_first=unlist(purrr::map2(.data$results_list,.data$first_result, ~match(.y,unlist(.x)))),
                        matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_first])))+1,)
      } else if(l=="last"){
        phenotype_extracted <- phenotype_extraction %>%
          dplyr::mutate(last_result= unlist(purrr::map(.data$results_list, ~ dplyr::last(stats::na.omit(.x)))),
                        position_last=unlist(purrr::map2(.data$results_list,.data$last_result, ~match(.y,unlist(.x)))),
                        matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_last])))+1,)
      } else if(l=="min"){
        phenotype_extracted <- phenotype_extraction %>%
          dplyr::mutate(min_result= unlist(purrr::map(.data$results_list, ~ min(unlist(.x),na.rm = T))),
                        position_min=unlist(purrr::map2(.data$results_list,.data$min_result, ~match(.y,unlist(.x)))),
                        matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_min])))+1,)
      } else if(l=="max"){
        phenotype_extracted <- phenotype_extraction %>%
          dplyr::mutate(max_result= unlist(purrr::map(.data$results_list, ~ max(unlist(.x),na.rm = T))),
                        position_max=unlist(purrr::map2(.data$results_list,.data$max_result, ~match(.y,unlist(.x)))),
                        matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_max])))+1,)
      }
      phenotype_extracted_complete <- phenotype_extracted %>%
        dplyr::select(.data$eid,tidyselect::all_of(position_col)) %>%
        dplyr::rename(eid=1,pheno=2,position_match=3) %>%
        tidyr::drop_na()
      # create either date and age column or just date column depending on input
      if(is.na(k)){
        phenotype_extraction_final <- phenotype_extracted_complete %>%
          dplyr::left_join(date_col, by="eid") %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
          dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                        earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
          dplyr::select(.data$eid,.data$pheno,.data$earliest_date) %>%
          dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
          purrr::set_names(c("eid",a,"earliest_date"))
        return(phenotype_extraction_final)
      } else {
        # extract age column
        age_col <- min_data %>%
          dplyr::select(.data$eid,tidyselect::matches(k))
        # vector used later
        age_col_names <- colnames(age_col)[-1]
        # add in age specific column
        phenotype_extraction_final <- phenotype_extracted_complete %>%
          dplyr::left_join(date_col, by="eid") %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
          dplyr::left_join(age_col, by="eid") %>%
          dplyr::mutate(date_list=purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                        earliest_date = purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
          dplyr::mutate(age_list=purrr::transpose(dplyr::across(tidyselect::all_of(age_col_names))),
                        pheno_age = unlist(purrr::map2(.data$age_list,.data$position_match, function(x,y) x[[y]]))) %>%
          dplyr::select(.data$eid,.data$pheno,.data$earliest_date,.data$pheno_age) %>%
          dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
          purrr::set_names(c("eid",a,"earliest_date",age_name))
      }
      return(phenotype_extraction_final)

    } else {
      # creating a variable used for naming a column
      PheWAS_ID_age <- paste0(a,"_age")
      # extract date column
      date_col <- min_data %>%
        dplyr::select(.data$eid,tidyselect::matches(j))
      # vector used later
      date_col_names <- colnames(date_col)[-1]
      # searching for all cols with the QC flag field ID
      flag_ID <- min_data %>%
        dplyr::select(.data$eid,tidyselect::matches(b))
      # can dplyr::select which QC flags to use as dplyr::filters if not discriminating then use 'all' in PheWAS_manifest
      if(n=="all") {
        just_QC <- flag_ID
      } else {
        # need to extract cases and controls separately, to achieve this need to identify unique values and set all non-cases or non-controls as NA
        all_values <- unique(as.vector(as.matrix(flag_ID[,-1])))
        flag_values <- as.numeric(unlist(strsplit(n,",")))
        # leave only values that represent non-accepted QC values
        na_values <- all_values[- which(all_values %in% flag_values)]
        # turn non-accepted qc flags into NAs
        just_QC <- per_data_field %>%
          dplyr::mutate(dplyr::across(dplyr::everything(),~replace(., . %in% na_values,NA)))
      }
      # converts any non-empty or NA string to the value of 1 as it is not important for this analysis why a flag is present.
      just_QC[,2:ncol(just_QC)][just_QC[,2:ncol(just_QC)] != ""] <- 1
      # converts to numeric
      just_QC <- just_QC %>%
        dplyr::mutate(dplyr::across(c(2:ncol(.)),as.numeric))

      # in some cases the QC flag has more measures than there are measures for the variable it is
      # flagging so there are more QC flags than measures, to account for this the flag measures
      # are combined into a single value per equivalent quantitative measure. They are then re-coded as the specific
      # cause of the QC flag is not relevant for our analysis.
      flag_cols_vector <- colnames(flag_ID)[-1]
      # use the vector of flag ID col names and identify how many visits this represents
      unique_visit <- paste0(b,unique(as.numeric(gsub(".*-(.+)\\..*", "\\1", flag_cols_vector))))
      # then dplyr::select each group of cols that represent each visit and convert into a single value if any value present
      split <- lapply(unique_visit, function (x) test <- just_QC %>%
                        dplyr::select(1,tidyselect::matches(x)))
      divide <- lapply(split, function (x) edit <- x %>%
                         dplyr::mutate(results = rowSums(dplyr::across(-.data$eid),na.rm = T)) %>%
                         dplyr::select(1,.data$results))
      divided <- divide %>%
        purrr::reduce(dplyr::left_join,by=c("eid")) %>%
        dplyr::mutate(filter = rowSums(dplyr::across(-.data$eid),na.rm = T)) %>%
        dplyr::filter(.data$filter>0) %>%
        dplyr::select(-.data$filter)
      # convert into a list that should be in equal length to the data_field the flag is used for
      flag_edit <- divided %>%
        dplyr::mutate(dplyr::across(tidyselect::all_of(colnames(.)[-1]), ~dplyr::na_if(.,0))) %>%
        dplyr::mutate(results_list=purrr::transpose(dplyr::across(tidyselect::all_of(colnames(.)[-1]))),
                      flags=purrr::map(.data$results_list, function(x) which(!is.na(unlist(x))))) %>%
        dplyr::select(.data$eid,.data$flags)
      # searching for all cols with the data_field
      per_data_field <- min_data %>%
        dplyr::select(.data$eid,tidyselect::matches(c))
      # vector needed later
      data_field_colnames <- colnames(per_data_field)[-1]
      # edits the original values based on the presence of a flag ID and then dplyr::selects the first valid value
      ID_edit <- per_data_field %>%
        tidyr::unite(col="results",2:length(colnames(.)),sep = ",",na.rm = T,remove = F) %>%
        dplyr::mutate(results=dplyr::na_if(.data$results,"")) %>%
        tidyr::drop_na(.data$results) %>%
        dplyr::left_join(flag_edit) %>%
        dplyr::mutate(results_list=purrr::transpose(dplyr::across(tidyselect::all_of(data_field_colnames))),
                      new_results=purrr::map2(.data$results_list,.data$flags, function(x,y) replace(unlist(x),list = unlist(y),values = NA)))
      # too computationally inefficient to create all variations so create only the one that is required
      if(l=="first"){
        phenotype_extracted <- ID_edit %>%
          dplyr::mutate(first_result= unlist(purrr::map(.data$new_results, ~ dplyr::first(stats::na.omit(.x)))),
                        position_first=unlist(purrr::map2(.data$new_results,.data$first_result, ~match(.y,unlist(.x)))),
                        matching_date_position_first=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_first])))+1,)
      } else if(l=="last"){
        phenotype_extracted <- ID_edit %>%
          dplyr::mutate(last_result= unlist(purrr::map(.data$new_results, ~ dplyr::last(stats::na.omit(.x)))),
                        position_last=unlist(purrr::map2(.data$new_results,.data$last_result, ~match(.y,unlist(.x)))),
                        matching_date_position_last=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_last])))+1,)
      } else if(l=="min"){
        phenotype_extracted <- ID_edit %>%
          dplyr::mutate(min_result= unlist(purrr::map(.data$new_results, ~ min(unlist(.x),na.rm = T))),
                        position_min=unlist(purrr::map2(.data$new_results,.data$min_result, ~match(.y,unlist(.x)))),
                        matching_date_position_min=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_min])))+1,)
      } else if(l=="max"){
        phenotype_extracted <- ID_edit %>%
          dplyr::mutate(max_result= unlist(purrr::map(.data$new_results, ~ max(unlist(.x),na.rm = T))),
                        position_max=unlist(purrr::map2(.data$new_results,.data$max_result, ~match(.y,unlist(.x)))),
                        matching_date_position_max=(as.numeric(gsub(".*-(.+)\\..*", "\\1", data_field_colnames[.data$position_max])))+1,)
      }
      phenotype_extracted_complete <- phenotype_extracted %>%
        dplyr::select(.data$eid,tidyselect::all_of(position_col)) %>%
        dplyr::rename(eid=1,pheno=2,position_match=3) %>%
        tidyr::drop_na()
      # used to add named_age or not as inputted
      if(is.na(k)){
        phenotype_extraction_final <- phenotype_extracted_complete %>%
          dplyr::left_join(date_col, by="eid") %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
          dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                        earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
          dplyr::select(.data$eid,.data$pheno,.data$earliest_date) %>%
          dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
          purrr::set_names(c("eid",a,"earliest_date"))
        return(phenotype_extraction_final)
      } else {
        # extract age column
        age_col <- min_data %>%
          dplyr::select(.data$eid,tidyselect::matches(k))
        # vector used later
        age_col_names <- colnames(age_col)[-1]
        # add in age specific column
        phenotype_extraction_final <- phenotype_extracted_complete %>%
          dplyr::left_join(date_col, by="eid") %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names), ~gsub(" .*", "", .x))) %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(date_col_names),as.character)) %>%
          dplyr::left_join(age_col, by="eid") %>%
          dplyr::mutate(date_list=purrr::transpose(dplyr::across(tidyselect::all_of(date_col_names))),
                        earliest_date = purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
          dplyr::mutate(age_list=purrr::transpose(dplyr::across(tidyselect::all_of(age_col_names))),
                        pheno_age = unlist(purrr::map2(.data$age_list,.data$position_match, function(x,y) x[[y]]))) %>%
          dplyr::select(.data$eid,.data$pheno,.data$earliest_date,.data$pheno_age) %>%
          dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
          purrr::set_names(c("eid",a,"earliest_date",age_name))
      }
      return(phenotype_extraction_final)
    }
  }
}

#' Makes combined data-field phenotypes from study data guided by PheWAS manifest
#'
#' @param a PheWAS_ID
#' @param b selecting whether first_last_max_min_value
#' @param c existence of limits or not
#' @param d lower_limit quant value
#' @param e upper_limit quant value
#' @param f data coding for age
#' @param PheWAS_manifest PheWAS_manifest
#' @param current_data_fields vector of current data fields available generated from previous code
#' @param data_field_variables dataframe of extracted data fields
#' @return A dataframe containing combined data-field phenotype
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data

combining_data_field <- function (a,b,c,d,e,f,PheWAS_manifest,current_data_fields,data_field_variables) {
  message(a)
  # selecting min or max values
  if(b=="min"){
    position_col <- c("min_result","position_min")
  } else if(b=="max"){
    position_col <- c("max_result","position_max")
  }
  # defining age column name is required
  age_name <- paste0(a,"_age")
  # value to search for
  phewas_ID_edit <- paste0(a,".")
  # combining values
  combinees_a <- PheWAS_manifest %>%
    dplyr::filter(stringr::str_detect(.data$PheWAS_ID,phewas_ID_edit)) %>%
    dplyr::pull(.data$PheWAS_ID)
  if(all(combinees_a %in% current_data_fields)){

    combining_variables <- data_field_variables[combinees_a] %>%
      purrr::reduce(dplyr::full_join,by="eid")
    #getting combinations of column names to split columns by
    ID_colnames <- colnames(combining_variables)[-1]
    ID_columns_no_date <- ID_colnames[!grepl("earliest_date", ID_colnames)]
    ID_colnames_just_date <- ID_colnames[!grepl(a, ID_colnames)]
    ID_columns_age <- ID_colnames[grepl("age", ID_colnames)]
    non_data_columns <- c(ID_colnames_just_date,ID_columns_age)
    # just the data
    ID_colnames_just_data <- ID_colnames [! ID_colnames %in% non_data_columns]
    # need to recode other columns that are not dates
    combined_variables <-combining_variables %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(ID_columns_no_date), as.double))
    # dplyr::selecting the component parts to recombine later
    just_data <- combined_variables %>%
      dplyr::select(.data$eid,tidyselect::all_of(ID_colnames_just_data))
    just_dates <- combined_variables %>%
      dplyr::select(.data$eid,tidyselect::all_of(ID_colnames_just_date)) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(ID_colnames_just_date), as.character))
    just_age <- combined_variables %>%
      dplyr::select(.data$eid,tidyselect::all_of(ID_columns_age))
    # creating columns to dplyr::select max and min and record position
    selecting_data <- just_data %>%
      dplyr::mutate(results_list=purrr::transpose(dplyr::across(tidyselect::all_of(ID_colnames_just_data))),
                    max_result=unlist(purrr::map(.data$results_list, ~ max(unlist(.x),na.rm = T))),
                    min_result=unlist(purrr::map(.data$results_list, ~ min(unlist(.x),na.rm = T))),
                    position_min=unlist(purrr::map2(.data$results_list,.data$min_result, ~match(.y,unlist(.x)))),
                    position_max=unlist(purrr::map2(.data$results_list,.data$max_result, ~match(.y,unlist(.x))))) %>%
      dplyr::select(.data$eid,tidyselect::all_of(position_col)) %>%
      dplyr::rename(eid=1,pheno=2,position_match=3)
    # option to add in age column as required, otherwise combine with just_dates and dplyr::select correct value
    if(f=="age"){
      phenotype_extraction_final <- selecting_data %>%
        dplyr::left_join(just_dates, by="eid") %>%
        dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(ID_colnames_just_date))),
                      earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
        dplyr::select(.data$eid,.data$pheno,.data$earliest_date) %>%
        dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
        purrr::set_names(c("eid",a,"earliest_date"))
      return(phenotype_extraction_final)
    } else if(f=="named"){
      phenotype_extraction_final <- selecting_data %>%
        dplyr::left_join(just_dates, by="eid") %>%
        dplyr::mutate(date_list = purrr::transpose(dplyr::across(tidyselect::all_of(ID_colnames_just_date))),
                      earliest_date=purrr::map2(.data$date_list,.data$position_match, function(x,y) x[[y]])) %>%
        dplyr::left_join(just_age, by="eid") %>%
        dplyr::mutate(age_list = purrr::transpose(dplyr::across(tidyselect::all_of(ID_columns_age))),
                      age=unlist(purrr::map2(.data$age_list,.data$position_match, function(x,y) x[[y]]))) %>%
        dplyr::select(.data$eid,.data$pheno,.data$earliest_date,.data$age) %>%
        dplyr::mutate(earliest_date=lubridate::ymd(.data$earliest_date)) %>%
        purrr::set_names(c("eid",a,"earliest_date",age_name))
      return(phenotype_extraction_final)
    }
  } else {return()}
}

#' Makes data-field phenotypes from study data guided by PheWAS manifest
#'
#' @param min_data Full path of the tab-separated file generated by the minimum_data.R script.
#' @param phenotype_save_file Full path for the save file for the generated concepts RDS to be used for phenotype creation.
#' @param N_cores Number of cores requested if parallel computing is desired. Defaults to single core computing.
#' @param PheWAS_manifest_overide Full path of an alternative to the in built file.
#' @param append_file Full path of an existing R object containing Data-Field phenotype information to which additional phenotypes derived by this script will be appended.
#' @return An RDS file containing lists of data-frames each representing a single data_field phenotype containing combined data-field phenotypes
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

data_field_phenotyping <- function(min_data,
                                   phenotype_save_file,
                                   N_cores,
                                   PheWAS_manifest_overide,
                                   append_file) {
# Load in defining variables
# min_data

if(!file.exists(min_data)){
  rlang::abort(paste0("'min_data' must be a file"))
}
min_data <- data.table::fread(min_data)

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

# Select if wanting to append existing data_field phenotype file
if(!is.null(append_file)){
  if(!file.exists(append_file)){
    rlang::abort(paste0("'append_file' must be a file"))
  }
  data_field_variables <- readRDS(append_file)
} else {
  data_field_variables <- list()
}
# save location
phenotype_save_location <- phenotype_save_file
new_folder <- stringr::str_remove(phenotype_save_file,basename(phenotype_save_file))
if(!dir.exists(new_folder)){
  dir.create(new_folder)}

# Data_field_phenotype function
# check to see if all the data is available and editing phenotype generation based on missing data and check to see if any existing phenotypes in data_field_variables
current_data_field <- gsub("-.*", "\\1", colnames(min_data))
current_IDs <- names(data_field_variables)
data_field <- PheWAS_manifest %>%
  tidyr::drop_na(.data$field_code) %>%
  dplyr::filter(.data$field_code %in% current_data_field,
                !.data$PheWAS_ID %in% current_IDs) %>%
  dplyr::mutate(search_id=paste0("^",.data$field_code,"-"),
                QC_flag_ID=ifelse(is.na(.data$QC_flag_ID),.data$QC_flag_ID,paste0("^",.data$QC_flag_ID,"-")),
                date_code=ifelse(is.na(.data$date_code),.data$date_code,paste0("^",.data$date_code,"-")),
                age_code=ifelse(is.na(.data$age_code),.data$age_code,paste0("^",.data$age_code,"-")),
                control_code=as.character(.data$control_code),
                case_code=as.character(.data$case_code))
# function to extract data_field phenotypes
if (is.numeric(N_cores)) {
  data_field_variables_created <- parallel::mcmapply(data_field_extraction,
                                                     data_field$PheWAS_ID,
                                                     data_field$QC_flag_ID,
                                                     data_field$search_id,
                                                     data_field$analysis,
                                                     data_field$case_code,
                                                     data_field$control_code,
                                                     data_field$limits,
                                                     data_field$lower_limit,
                                                     data_field$upper_limit,
                                                     data_field$date_code,
                                                     data_field$age_code,
                                                     data_field$first_last_max_min_value,
                                                     data_field$exclude_case_control_crossover,
                                                     data_field$QC_filter_values,
                                                     MoreArgs=list(min_data=min_data),
                                                     SIMPLIFY = F,
                                                     mc.cores = N_cores,
                                                     USE.NAMES = T)
} else {
  data_field_variables_created <- mapply(data_field_extraction,
                                         data_field$PheWAS_ID,
                                         data_field$QC_flag_ID,
                                         data_field$search_id,
                                         data_field$analysis,
                                         data_field$case_code,
                                         data_field$control_code,
                                         data_field$limits,
                                         data_field$lower_limit,
                                         data_field$upper_limit,
                                         data_field$date_code,
                                         data_field$age_code,
                                         data_field$first_last_max_min_value,
                                         data_field$exclude_case_control_crossover,
                                         data_field$QC_filter_values,
                                         MoreArgs=list(min_data=min_data),
                                         SIMPLIFY = F,
                                         USE.NAMES = T)
}
# append to existing list either empty or with previously created phenotypes
data_field_variables <- append(data_field_variables,data_field_variables_created)
data_field_variables <- purrr::compact(data_field_variables)
current_data_fields <- names(data_field_variables)

# Combined_data_field_phenotype function
data_field_comb <- PheWAS_manifest %>%
  dplyr::filter(.data$quant_combination==1)

combine_search <- paste0(data_field_comb$PheWAS_ID,".")

combinees <- PheWAS_manifest %>%
  dplyr::filter(stringr::str_detect(.data$PheWAS_ID,paste(combine_search,collapse = "|"))) %>%
  dplyr::pull(.data$PheWAS_ID)

if(isFALSE(any(combinees %in% current_data_fields))){

} else {
  # function
  if (length(data_field_comb >0)) {
    filed_ID_combined <- mapply(combining_data_field,
                                data_field_comb$PheWAS_ID,
                                data_field_comb$first_last_max_min_value,
                                data_field_comb$limits,
                                data_field_comb$lower_limit,
                                data_field_comb$upper_limit,
                                data_field_comb$age_col,
                                MoreArgs = list(PheWAS_manifest=PheWAS_manifest,
                                                current_data_fields=current_data_fields,
                                                data_field_variables=data_field_variables),
                                SIMPLIFY = F,
                                USE.NAMES = T)
    data_field_variables <- append(data_field_variables,filed_ID_combined)
    data_field_variables <- purrr::compact(data_field_variables)
  }
}
# save file
saveRDS(data_field_variables,phenotype_save_location)

}
