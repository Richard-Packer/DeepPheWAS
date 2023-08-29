

#' Searches and extracts named columns from data file
#'
#' @param x full file path to data file
#' @param y vector of column names to extract
#' @param z dataframe of exclusions single column no header
#' @return A dataframe containing only columns of interest
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data

tab_data_clean <- function (x,y,z) {
  # finding column names in data set

  # need a warning message here that

  tabdata_colname <- colnames(data.table::fread(x, nrows=0))
  colname_search <- stringr::str_subset(tabdata_colname,paste0(y,collapse = "|"))
  colname_search <- colname_search[!is.na(colname_search)]
  # and extracting only required columns
  tabdata <- data.table::fread(x,select = colname_search) %>%
    dplyr::rename(eid=1) %>%
    dplyr::filter(!.data$eid %in% z$V1)
  return(tabdata)
}

#' Searches and extracts named columns from data file
#'
#' @param data_folder Full path of the directory that contains the data files that will be formatted and concatenated.
#' @param data_files Comma separated full file paths of the data files that will be formatted and concatenated.
#' @param r_format Specific to UK Biobank data. Specify if the input for the data have been downloaded using the R option. If the data have been downloaded using the .csv or .txt options, then no input is required (default).
#' @param data_field_ID Full path to the file containing the field_IDs required for Deep-PheWAS if not using default file. Is a plain text file with no header one field-ID per row. Field-ID is a numeric value, example field-ID 54 is UK biobank assessment centre.
#' @param data_name_pattern Character string for isolating data files if these are in a directory with other files. Defaults to using all files in the directory specified by the data_folder argument.
#' @param N_cores Number of cores requested if parallel computing is desired. Defaults to single core computing.
#' @param save_loc Full path to save file location for minimum_tab_data.
#' @param exclusions Full path to the file containing individuals to be excluded from the analysis. Defaults behaviour is to retain all individuals.
#' @return A dataframe containing only columns of interest
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

minimum_data_R <- function(data_folder,
                           data_files,
                           r_format,
                           data_field_ID,
                           data_name_pattern,
                           N_cores,
                           save_loc,
                           exclusions){


  # Evaluating arguments ----------------------------------------------------
  # data_folder
  if(!is.null(data_folder)){
    if(!dir.exists(data_folder)){
      rlang::abort(paste0("'data_folder' must be an active directory"))
    }
  }
  # data_files
  if(!is.null(data_files)){
    tab_files <- unlist(stringr::str_split(data_files,pattern = ","))
    checking_files <- unlist(lapply(tab_files,function(x) file.exists(x)))
    if(!all(checking_files)){
      not_files <- tab_files[which(checking_files==F)]
      rlang::abort(paste0("all 'data_files' must be files ",not_files," does not appear to be a readable file"))
    }
  }
  #exclusions
  if(!is.null(exclusions)){
    if(!file.exists(exclusions))
    rlang::abort(paste0("'exclusions' must be a readable file"))
  }
  # save_loc
  if(!dir.exists(dirname(save_loc))){
    dir.create(dirname(save_loc),recursive = T)
  }

  # Load in and define variables -----------------------------------------------------------------
  # N_cores
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
  # required data_fields
  if(data_field_ID!="fields-minimum.txt") {
    field_ID_headings <- data.table::fread(data_field_ID, header = F) } else {
      field_ID_headings <- data.table::fread(system.file("extdata","fields-minimum.txt.gz", package = "DeepPheWAS"), header = F)
    }
  adding_eid <- data.frame(conversion=c("eid","f.eid"))
  # selecting data format for UK-Biobank only
  if(r_format){
    headings <- field_ID_headings %>%
      dplyr::mutate(conversion=paste0("f.",.data$V1,".")) %>%
      dplyr::bind_rows(adding_eid) %>%
      dplyr::pull(.data$conversion)
  } else {
    headings <- field_ID_headings %>%
      dplyr::mutate(conversion=paste0("^",.data$V1,"-")) %>%
      dplyr::bind_rows(adding_eid) %>%
      dplyr::pull(.data$conversion)
  }
  # exclusions
  if(!is.null(exclusions)) {
    exclusions <- data.table::fread(exclusions, header = F)
  } else {
    exclusions <- data.frame(V1=NA)
  }
  # selecting files to search
  if(!is.null(data_folder)) {
    if(is.null(data_name_pattern)) {
      tab_files <- list.files(path = data_folder, full.names = T)
      checking_files <- unlist(lapply(tab_files,function(x) file.exists(x)))
      if(!all(checking_files)){
        not_files <- tab_files[which(checking_files==F)]
        rlang::abort(paste0(not_files," in ",data_folder," are not readable all files in ",data_folder," must be readable"))
      }
    } else {
      tab_files <- list.files(path = data_folder, pattern = data_name_pattern, full.names = T)
      checking_files <- unlist(lapply(tab_files,function(x) file.exists(x)))
      if(!all(checking_files)){
        not_files <- tab_files[which(checking_files==F)]
        rlang::abort(paste0(not_files," in ",data_folder," with pattern ",data_name_pattern," are not readable all files in ",data_folder," with pattern ",data_name_pattern," must be readable"))
      }
    }
  } else if(!is.null(data_files)){
    tab_files <- unlist(stringr::str_split(data_files,pattern = ","))
  }
  # Run functions -------------------------------------------------------
  if(is.numeric(N_cores)){
    min_tab_data <- parallel::mcmapply(tab_data_clean,
                                       tab_files,
                                       MoreArgs = list(y=headings, z=exclusions),
                                       mc.cores = N_cores,
                                       SIMPLIFY = F)
  } else {
    min_tab_data <- mapply(tab_data_clean,
                           tab_files,
                           MoreArgs = list(y=headings, z=exclusions),
                           SIMPLIFY = F)
  }
  # Use purr to map and join Worst case scenario numbers don't match but it just introduces NA's.
  minimum_tab <- min_tab_data %>%
    purrr::reduce(dplyr::full_join, by="eid")
  if(r_format) {
    colnames(minimum_tab) <- stringr::str_replace(colnames(minimum_tab),"f\\.","")
    colnames(minimum_tab) <- stringr::str_replace(colnames(minimum_tab),"\\.","-")
  }
  data.table::fwrite(minimum_tab, save_loc, sep = "\t", na = NA)
}
