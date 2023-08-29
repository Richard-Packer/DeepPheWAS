#' Simple function that extracts names of lists within lists.
#'
#' @param x Name of the first level of list.
#' @param y list where names are extracted.
#' @return a list of names from each internal list
#' @keywords internal
inner_names <- function(x,y){
  all_ID <- y[[x]] %>%
    dplyr::pull(.data$PheWAS_ID)
}


#' Simple function that adds to existing list names of lists within lists.
#'
#' @param a Name of list to append.
#' @param b RDS of original results.
#' @param c new results
#' @param d name of the trait
#' @param e group name
#' @return an appended list.
#' @keywords internal
adding_to_results <- function(a,b,c,d,e){
  original_table <- b[[a]]
  new_results <- c[[a]]

  combined_results <- original_table %>%
    dplyr::bind_rows(new_results) %>%
    dplyr::mutate("name"={{d}},
                  table_save=.data$name,
                  name_group=paste0(.data$name,"_",{{e}}))
}

#' Applies association analysis per grouping variable for non-GRS data.
#'
#' @param x Group name
#' @param y Full file location of a phenotype file
#' @param tested_variables genetic variable data to be analysed from the non_GRS data.
#' @param GRS_input Data frame guide to location of GRS data, and phenotype data per trait/grouping combination.
#' @param non_GRS_data Data frame of non-GRS data, eid column, with each column representing a genetic variable.
#' @param covariates Full file path to the covariate file edited for use in association testing, see user guide.
#' @param all_phenos Data frame containing all possible phenotypes including age column covariates.
#' @param PheWAS_manifest Data frame of the PheWAS_manifest containing information on phenotypes.
#' @param analysis_save_name Core save name for the analysis.
#' @param quant_min_cases minimum number of quantitative cases.
#' @param binary_min_cases minimum number of binary cases.
#' @param covariate_col_names column names of covariates.
#' @param no_covariate T or F if no covariate file.
#' @return non-GRS results as list item.
#' @keywords internal
per_group_none_GRS  <- function(x,y,
                                tested_variables,
                                non_GRS_data,
                                GRS_input,
                                covariates,
                                all_phenos,
                                PheWAS_manifest,
                                analysis_save_name,
                                quant_min_cases,
                                binary_min_cases,
                                covariate_col_names,
                                no_covariate){

  GRS_results <- mapply(GRS_association,
                        a=tested_variables,
                        d=tested_variables,
                        MoreArgs = list(b=non_GRS_data,c=y,e=x,
                                        GRS_input=GRS_input,
                                        non_GRS_data=non_GRS_data,
                                        covariates=covariates,
                                        all_phenos=all_phenos,
                                        PheWAS_manifest=PheWAS_manifest,
                                        analysis_save_name=analysis_save_name,
                                        quant_min_cases=quant_min_cases,
                                        binary_min_cases=binary_min_cases,
                                        covariate_col_names=covariate_col_names,
                                        no_covariate=no_covariate),
                        SIMPLIFY = F,
                        USE.NAMES = T) %>%
    purrr::reduce(dplyr::bind_rows)
  return(GRS_results)
}

#' Applies association analysis per grouping variable for GRS data.
#'
#' @param x Group name
#' @param GRS_input Data frame guide to location of GRS data, and phenotype data per trait/grouping combination.
#' @param non_GRS_data Data frame of non-GRS data, eid column, with each column representing a genetic variable.
#' @param covariates Full file path to the covariate file edited for use in association testing, see user guide.
#' @param all_phenos Data frame containing all possible phenotypes including age coloumn covariates.
#' @param PheWAS_manifest Data frame of the PheWAS_manifest containing information on phenotypes.
#' @param analysis_save_name Core save name for the analysis.
#' @param quant_min_cases minimum number of quantitative cases.
#' @param binary_min_cases minimum number of binary cases.
#' @param covariate_col_names column names of covariates.
#' @param old_results existing results as list item.
#' @param no_covariate T or F if no covariate file.
#' @return GRS results as a saved RDS object.
#' @keywords internal
per_group_per_trait_GRS <- function(x,
                                    GRS_input,
                                    non_GRS_data,
                                    covariates,
                                    all_phenos,
                                    PheWAS_manifest,
                                    analysis_save_name,
                                    quant_min_cases,
                                    binary_min_cases,
                                    covariate_col_names,
                                    old_results,
                                    no_covariate) {

  selected_trait <- x

  save_name <- unique(selected_trait$save_location)

  GRS_results <- mapply(GRS_association,
                        selected_trait$group,
                        selected_trait$trait,
                        selected_trait$genetic_data,
                        selected_trait$phenotype_data,
                        selected_trait$column_name,
                        MoreArgs = list(GRS_input=GRS_input,
                                        non_GRS_data=non_GRS_data,
                                        covariates=covariates,
                                        all_phenos=all_phenos,
                                        PheWAS_manifest=PheWAS_manifest,
                                        analysis_save_name=analysis_save_name,
                                        quant_min_cases=quant_min_cases,
                                        binary_min_cases=binary_min_cases,
                                        covariate_col_names=covariate_col_names,
                                        no_covariate=no_covariate),
                        SIMPLIFY = F,
                        USE.NAMES = T)


  if(is.null(old_results)){
    saveRDS(GRS_results,save_name)
  } else {
    to_change <- old_results[names(GRS_results)]

    updated_results <- mapply(adding_to_results,names(to_change),MoreArgs = list(b=old_results,c=GRS_results,d=selected_trait$trait,e=selected_trait$group), SIMPLIFY = F)

    unchanged <- old_results[names(old_results)[!names(old_results) %in% names(updated_results)]]
    if(length(unchanged)>0){
    final <- append(unchanged,updated_results)[names(old_results)]
    } else {
  final <- updated_results
    }
    saveRDS(final,save_name)
  }
}

#' Simple helper function for making folders
#'
#' @param x Per trait map of files
#' @param analysis_folder Location where analysis output is saved to.
#' @return Plink association results.
#' @keywords internal
making_GRS_folders <- function(x,analysis_folder) {

  GRS_per_trait <- x
  trait_name <- unique(x$trait)

  lapply(GRS_per_trait$group, function(x) dir.create(paste0(analysis_folder,"/",trait_name,"/",x),recursive = T))

}

#' Runs R association testing.
#'
#' @param a PheWAS_ID
#' @param b age column information whether generic age covariate or age at measurement extracted from the phenotype files.
#' @param c analysis option binary or quantitiative
#' @param d variable name
#' @param phenotypes data frame of all phenotypes for current grouping.
#' @param GRS_cov covariates for association testing
#' @param PheWAS_manifest Data frame of the PheWAS_manifest containing information on phenotypes.
#' @param quant_min_cases minimum cases for quantitative phenotypes
#' @param binary_min_cases minimum cases for binary phenotypes
#' @param covariate_col_names column names of covariates.
#' @param all_phenos Data frame containing all possible phenotypes including age coloumn covariates.
#' @return A dataframe of association results per phenotype.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
GRS_association_testing <- function (a,b,c,d,phenotypes,PheWAS_manifest,GRS_cov,quant_min_cases,binary_min_cases,covariate_col_names,all_phenos) {
  extracted_cols <- c(a,b)
message(paste0(a))
  phenotype_extracted <- phenotypes %>%
    dplyr::select(.data$eid,tidyselect::any_of(extracted_cols)) %>%
    dplyr::left_join(GRS_cov) %>%
    tidyr::drop_na()

  phenotype_selected <- phenotype_extracted %>%
    dplyr::select({{a}}) %>%
    dplyr::rename(selected_phenotype=1)

  type <- all_phenos %>%
    dplyr::filter(.data$PheWAS_ID=={{a}}) %>%
    dplyr::select(.data$analysis) %>%
    dplyr::pull()

  if (type=="quant"){
    if(nrow(phenotype_selected)<quant_min_cases){
      return()
    }

  } else if(type=="binary") {
    if(sum(phenotype_selected$selected_phenotype)<binary_min_cases)
      return()
  }

  if(b=="") {
    formula_covariates <- paste(colnames(covariate_col_names), collapse = " + ")
  } else {
    formula_covariates <- paste0(" + ", paste(c(b,colnames(covariate_col_names)), collapse = " + "))
  }

  if (c=="binomial") {
    my_form <- paste0("as.factor(",a,")~ ",d,formula_covariates)
    my_formula <- stats::as.formula(my_form)

    my_glm <- stats::glm(my_formula, family = c, data = phenotype_extracted, na.action = stats::na.omit)

    p_value <- (stats::coef(summary(my_glm))[,4])[d]
    SE <- unname((summary(my_glm)$coefficients[, 2])[d])
    P <- unname(p_value)
    if (is.na(P)) {
      P <- NA
      odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>%
        dplyr::rename(OR=1,L95=2,U95=3,SE=4)

    } else if (P > 0.05) {
      OR_raw <- exp(stats::coef(my_glm))[d]
      OR <- unname(OR_raw)
      if (is.na(OR)) {
        odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>%
          dplyr::rename(OR=1,L95=2,U95=3,SE=NA)
      } else {
        odds <- exp(stats::coef(my_glm)[d])
        odds_df <- data.frame(list(OR=unname(odds),L95=NA,U95=NA,SE=SE)) %>%
          dplyr::rename(OR=1,L95=2,U95=3,SE=4)
      }
    } else {
      OR_raw <- exp(stats::coef(my_glm))[d]
      OR <- unname(OR_raw)
      if (is.na(OR)) {
        odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>%
          dplyr::rename(OR=1,L95=2,U95=3,SE=4)
      } else {
        odds_CI <- exp(odds_CI <- c(Beta = stats::coef(my_glm)[d], stats::confint(my_glm,d),SE=SE))
        odds_df <- data.frame(as.list(odds_CI)) %>%
          dplyr::rename(OR=1,L95=2,U95=3,SE=4)
      }
    }


  } else if (c=="gaussian") {
    my_form <- paste0(a," ~ ",d,formula_covariates)
    my_formula <- stats::as.formula(my_form)

    my_glm <- stats::glm(my_formula, family = c, data = phenotype_extracted, na.action = stats::na.omit)

    p_value <- (stats::coef(summary(my_glm))[,4])[d]
    SE <- unname((summary(my_glm)$coefficients[, 2])[d])
    P <- unname(p_value)
    if (is.na(P)) {
      P <- NA
      odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>%
        dplyr::rename(Beta=1,L95=2,U95=3,SE=4)

    } else if (P > 0.05) {
      OR_raw <- stats::coef(my_glm)[d]
      OR <- unname(OR_raw)
      if (is.na(OR)) {
        odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>%
          dplyr::rename(Beta=1,L95=2,U95=3,SE=4)
      } else {
        odds <- stats::coef(my_glm)[d]
        odds_df <- data.frame(list(Beta=unname(odds),L95=NA,U95=NA,SE=SE)) %>%
          dplyr::rename(Beta=1,L95=2,U95=3,SE=4)
      }
    } else {
      OR_raw <- stats::coef(my_glm)[d]
      OR <- unname(OR_raw)
      if (is.na(OR)) {
        odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>%
          dplyr::rename(Beta=1,L95=2,U95=3,SE=4)
      } else {
        odds_CI <- c(Beta = stats::coef(my_glm)[d], stats::confint(my_glm,d),SE=SE)
        odds_df <- data.frame(as.list(odds_CI)) %>%
          dplyr::rename(Beta=1,L95=2,U95=3,SE=4)
      }
    }
  }

  P_df <- data.frame(P)

  name_df <- data.frame(a) %>%
    dplyr::rename(PheWAS_ID=1)

  df_results <- name_df %>%
    dplyr::bind_cols(P_df,odds_df)

  return(df_results)

}
#' Runs R association testing per group per trait.
#'
#' @param e group identifier.
#' @param a Trait name.
#' @param b Full file path to location of the genetic data.
#' @param c Full file path to location of the phenotype data.
#' @param d column name within the genetic data to be extracted for association analysis.
#' @param GRS_input Data frame guide to location of GRS data, and phenotype data per trait/grouping combination.
#' @param non_GRS_data Data frame of non-GRS data, eid column, with each column representing a genetic variable.
#' @param covariates Data frame of covariates for analysis.
#' @param all_phenos Data frame containing all possible phenotypes including age coloumn covariates.
#' @param PheWAS_manifest Data frame of the PheWAS_manifest containing information on phenotypes.
#' @param analysis_save_name Core save name for the analysis.
#' @param quant_min_cases minimum number of quantitative cases.
#' @param binary_min_cases minimum number of binary cases.
#' @param covariate_col_names column names of covariates.
#' @param no_covariate T or F if no covariate file.
#' @return Plink association results.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
GRS_association <- function(e,a,b,c,d,
                            GRS_input,
                            non_GRS_data,
                            covariates,
                            all_phenos,
                            PheWAS_manifest,
                            analysis_save_name,
                            quant_min_cases,
                            binary_min_cases,
                            covariate_col_names,
                            no_covariate) {

  results_name <- a

  if(!is.null(GRS_input)){
    table_save_name <- results_name
  } else if (!is.null(non_GRS_data)){
    table_save_name <- analysis_save_name
  }

  variable_name <- d

  ## Load the GRS

  GRS <- data.table::fread(b) %>%
    dplyr::rename(eid=1)
  if(!is.null(GRS_input)){
  # combine GRS with covariates
  GRS_cov <- GRS %>%
    dplyr::left_join(covariates) %>%
    tidyr::drop_na()
} else if(!is.null(non_GRS_data)){
  GRS_cov <- GRS %>%
    dplyr::left_join(covariates) %>%
    dplyr::select("eid",tidyr::any_of(a)) %>%
    tidyr::drop_na()
  }
  # link phenotypes with GRS
  phenotypes <- data.table::fread(c) %>%
    dplyr::select(1,tidyselect::any_of(all_phenos$PheWAS_ID))
  available_phenotypes <- colnames(phenotypes)

  if(no_covariate) {
    association_guide <- all_phenos %>%
      dplyr::filter(.data$PheWAS_ID %in% available_phenotypes) %>%
      tidyr::drop_na(.data$analysis) %>%
      dplyr::mutate(age_col_complete = "",
                    analysis_option = ifelse(.data$analysis=="quant","gaussian","binomial")) %>%
      dplyr::select(.data$PheWAS_ID,.data$age_col_complete,.data$analysis_option)
  } else {
    association_guide <- all_phenos %>%
      dplyr::filter(.data$PheWAS_ID %in% available_phenotypes) %>%
      tidyr::drop_na(.data$analysis) %>%
      dplyr::mutate(temp_name=stringr::str_remove(.data$PheWAS_ID,"_male|_female|_age_of_onset"),
                    age_col_complete = ifelse(.data$age_col=="named",paste0(.data$temp_name,"_age"),"age"),
                    analysis_option = ifelse(.data$analysis=="quant","gaussian","binomial")) %>%
      dplyr::select(.data$PheWAS_ID,.data$age_col_complete,.data$analysis_option)
  }


  ## all association results have tables with Phewas_ID,P,OR,Beta
  association_results <- mapply(GRS_association_testing,
                                association_guide$PheWAS_ID,
                                association_guide$age_col_complete,
                                association_guide$analysis_option,
                                MoreArgs = list(d=variable_name,
                                                phenotypes=phenotypes,
                                                PheWAS_manifest=PheWAS_manifest,
                                                GRS_cov=GRS_cov,
                                                quant_min_cases=quant_min_cases,
                                                binary_min_cases=binary_min_cases,
                                                covariate_col_names=covariate_col_names,
                                                all_phenos=all_phenos),
                                SIMPLIFY = F) %>%
    purrr::reduce(dplyr::bind_rows)

  if ("Beta" %in% colnames(association_results) & "OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>%
      dplyr::rowwise() %>%
      dplyr::mutate(P=ifelse(.data$P>0.05,signif(.data$P, 2),signif(.data$P, 3)),
                    L95=ifelse(is.na(.data$Beta),signif(.data$L95, 5),signif(.data$L95, 2)),
                    U95=ifelse(is.na(.data$Beta),signif(.data$U95, 5),signif(.data$U95, 2)),
                    Beta=signif(.data$Beta, 2),
                    name=results_name,
                    SE=signif(.data$SE, 2),
                    effect_direction=ifelse(is.na(.data$OR),ifelse(is.na(.data$Beta),NA,ifelse(.data$Beta>0,"positive","negative")),ifelse(.data$OR>1,"positive","negative")),
                    OR=signif(.data$OR, 5),
                    group=e,
                    name_group=paste0(.data$name,"_",.data$group),
                    table_save=table_save_name,
                    join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female|_age_of_onset"),
                    remove_string=paste0(.data$join_name,"_"),
                    sex_pheno=stringr::str_remove(.data$PheWAS_ID,.data$remove_string)) %>%
      dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
      dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc),
                    phenotype=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),paste0("Age of onset of ",.data$phenotype),
                                     ifelse(stringr::str_detect(.data$PheWAS_ID,"_male"),paste0("Male stratified variant of ",.data$phenotype),
                                            ifelse(stringr::str_detect(.data$PheWAS_ID,"_female"),paste0("Female stratified variant of ",.data$phenotype),
                                                   .data$phenotype)))) %>%
  dplyr::select(.data$name,.data$PheWAS_ID,.data$phenotype,.data$P,.data$OR,.data$Beta,.data$L95,.data$U95,.data$SE,phenotype_group=.data$pheno_group,.data$group_narrow,.data$short_desc,.data$effect_direction,.data$group,.data$name_group,.data$table_save)

  } else if ("Beta" %in% colnames(association_results) & !"OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>%
      dplyr::rowwise() %>%
      dplyr::mutate(P=ifelse(.data$P>0.05,signif(.data$P, 2),signif(.data$P, 3)),
                    L95=signif(.data$L95, 2),
                    U95=signif(.data$U95, 2),
                    Beta=signif(.data$Beta, 2),
                    name=results_name,
                    SE=signif(.data$SE, 2),
                    effect_direction=ifelse(is.na(.data$Beta),NA,ifelse(.data$Beta>0,"positive","negative")),
                    group=e,
                    name_group=paste0(.data$name,"_",.data$group),
                    table_save=table_save_name,
                    join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female|_age_of_onset"),
                    remove_string=paste0(.data$join_name,"_"),
                    sex_pheno=stringr::str_remove(.data$PheWAS_ID,.data$remove_string)) %>%
      dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
      dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc),
                    phenotype=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),paste0("Age of onset of ",.data$phenotype),
                                     ifelse(stringr::str_detect(.data$PheWAS_ID,"_male"),paste0("Male stratified variant of ",.data$phenotype),
                                            ifelse(stringr::str_detect(.data$PheWAS_ID,"_female"),paste0("Female stratified variant of ",.data$phenotype),
                                                   .data$phenotype)))) %>%
      dplyr::select(.data$name,.data$PheWAS_ID,.data$phenotype,.data$P,.data$Beta,.data$L95,.data$U95,.data$SE,phenotype_group=.data$pheno_group,.data$group_narrow,.data$short_desc,.data$effect_direction,.data$group,.data$name_group,.data$table_save)

  } else if (!"Beta" %in% colnames(association_results) & "OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>%
      dplyr::rowwise() %>%
      dplyr::mutate(P=ifelse(.data$P>0.05,signif(.data$P, 2),signif(.data$P, 3)),
                    L95=signif(.data$L95, 5),
                    U95=signif(.data$U95, 5),
                    name=results_name,
                    SE=signif(.data$SE, 2),
                    effect_direction=ifelse(is.na(.data$OR),NA,ifelse(.data$OR>1,"positive","negative")),
                    OR=signif(.data$OR, 5),
                    group=e,
                    name_group=paste0(.data$name,"_",.data$group),
                    table_save=table_save_name,
                    join_name = stringr::str_remove(.data$PheWAS_ID,"_male|_female|_age_of_onset"),
                    remove_string=paste0(.data$join_name,"_"),
                    sex_pheno=stringr::str_remove(.data$PheWAS_ID,.data$remove_string)) %>%
      dplyr::left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>%
      dplyr::mutate(short_desc=ifelse(.data$sex_pheno!=.data$PheWAS_ID,paste0(.data$short_desc," (",.data$sex_pheno,")"), .data$short_desc),
                    phenotype=ifelse(stringr::str_detect(.data$PheWAS_ID,"_age_of_onset"),paste0("Age of onset of ",.data$phenotype),
                                     ifelse(stringr::str_detect(.data$PheWAS_ID,"_male"),paste0("Male stratified variant of ",.data$phenotype),
                                            ifelse(stringr::str_detect(.data$PheWAS_ID,"_female"),paste0("Female stratified variant of ",.data$phenotype),
                                                   .data$phenotype)))) %>%
      dplyr::select(.data$name,.data$PheWAS_ID,.data$phenotype,.data$P,.data$OR,.data$L95,.data$U95,.data$SE,phenotype_group=.data$pheno_group,.data$group_narrow,.data$short_desc,.data$effect_direction,.data$group,.data$name_group,.data$table_save)

  }

  return(association_results_edit)

}

#' Performs association testing in R for either genetic risk score (GRS) and non-GRS data
#'
#' @param analysis_folder Full path to the folder location of the folder hosting this analysis.
#' @param phenotype_files Comma separated full file paths of the phenotype data. For example /home/phenotypes/EUR_pheno,/home/phenotypes/AFR_pheno
#' @param covariates Full file path to the covariate file edited for use in association testing, see user guide.
#' @param GRS_input Full file path of the GRS_input csv file. See user guide for more information on format.
#' @param non_GRS_data Full file path to the non_GRS genetic data. See user guide for more information on format.
#' @param group_name_overide Comma-separated list containing alternative group names. By default, group names are extracted from the suffix of the file names provided in the phenotype_files argument. For example, if a file named /home/phenotypes/EUR_phenotypes.csv were provided, the corresponding group name would be "EUR". This argument allows for a different group name to be specified. The order of names provided should match the files specified in the phenotype_files argument.
#' @param PheWAS_manifest_overide Full file path of the alternative PheWAS_manifest file.
#' @param analysis_name Name for the analysis, is used later in saving tables, so should distinguish between other analyses. For GRS analysis this name is always the trait being analysed.
#' @param N_cores Number of cores requested if wanting to use parallel computing.
#' @param phenotype_inclusion_file Full file path to a plain txt file containing single column NO header containing full PheWAS_ID of phenotypes that will be included. Cannot be used with phenotype_exclusion_file argument.
#' @param phenotype_exclusion_file Full file path to a plain txt file containing single column NO header containing full PheWAS_ID of phenotypes that will be excluded. Cannot be used with phenotype_exclusion_file argument.
#' @param binary_Case_N Number that represents the minimum number of cases for binary phenotype inclusion. Default=50
#' @param quantitative_Case_N Number that represents the minimum number of cases for quantitative phenotype inclusion. Default=100
#' @param use_existing_results_file Full file path of existing results. Will append that file and save new file in save location inputted in the rest of commands. Speeds up time of analysis by utilising existing results rather than re-running all associations.
#' @return Association test results as a RDS object.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
R_association_testing <- function(analysis_folder,
                                  phenotype_files,
                                  covariates,
                                  GRS_input,
                                  non_GRS_data,
                                  group_name_overide,
                                  PheWAS_manifest_overide,
                                  analysis_name,
                                  N_cores,
                                  phenotype_inclusion_file,
                                  phenotype_exclusion_file,
                                  binary_Case_N,
                                  quantitative_Case_N,
                                  use_existing_results_file){
. <- NULL
  # read in files
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
  if(is.na(as.numeric(quantitative_Case_N))){
    rlang::abort(paste0("'quantitative_Case_N' must be a numeral"))
  }
  quant_min_cases <- as.numeric(quantitative_Case_N)
  if(is.na(as.numeric(binary_Case_N))){
    rlang::abort(paste0("'binary_Case_N' must be a numeral"))
  }
  binary_min_cases <-as.numeric(binary_Case_N)

  if(!is.null(covariates)) {
    if(!file.exists(covariates)){
      rlang::abort(paste0("'covariates' must be a file"))
    }
    covariates <- data.table::fread(covariates)
    covariate_col_names <- covariates %>%
      dplyr::select(-.data$eid,-.data$age)
    no_covariate <- FALSE
  } else {
    covariates <- data.frame(eid=NA)
    covariate_col_names <- covariates %>%
      dplyr::select(-.data$eid)
    no_covariate <- TRUE
  }

  # defining all_phenos
  additional_phenos <- PheWAS_manifest %>%
    dplyr::filter(.data$included_in_analysis==1) %>%
    dplyr::mutate(male=paste0(.data$PheWAS_ID,"_male"),
                  female=paste0(.data$PheWAS_ID,"_female"),
                  age_of_onset=paste0(.data$PheWAS_ID,"_age_of_onset"),
                  male_analysis=.data$analysis,
                  female_analysis=.data$analysis,
                  age_of_onset_analysis="quant",
                  male_age_col=.data$age_col,
                  female_age_col=.data$age_col,
                  age_of_onset_age_col="age") %>%
    dplyr::select(.data$PheWAS_ID,.data$male,.data$female,.data$analysis,.data$male_analysis,.data$female_analysis,.data$age_col,.data$male_age_col,.data$female_age_col,.data$age_of_onset,.data$age_of_onset_analysis,.data$age_of_onset_age_col) %>%
    dplyr::mutate(age_col_complete = ifelse(.data$age_col=="named",paste0(.data$PheWAS_ID,"_age"),"age"))

  all_possible_phenotypes <- data.frame(PheWAS_ID=c(additional_phenos$PheWAS_ID,additional_phenos$male,additional_phenos$female,additional_phenos$age_of_onset),
                                        analysis=c(additional_phenos$analysis,additional_phenos$male_analysis,additional_phenos$female_analysis,additional_phenos$age_of_onset_analysis),
                                        age_col=c(additional_phenos$age_col,additional_phenos$male_age_col,additional_phenos$female_age_col,additional_phenos$age_of_onset_age_col))

  if(!is.null(phenotype_inclusion_file)){
    if(!file.exists(phenotype_inclusion_file)){
      rlang::abort(paste0("'phenotype_inclusion_file' must be a file"))
    }
    phenotype_inclusion_ID <- data.table::fread(phenotype_inclusion_file, header = F) %>%
      dplyr::pull(1)
    nearly_all_phenos <- all_possible_phenotypes %>%
      dplyr::filter(.data$PheWAS_ID %in% phenotype_inclusion_ID)

    age_phenos <- nearly_all_phenos %>%
      dplyr::mutate(ID2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female|_age_of_onset"),NA,.data$PheWAS_ID)) %>%
      tidyr::drop_na(.data$ID2) %>%
      dplyr::mutate(age_col_complete = ifelse(.data$age_col=="named",paste0(.data$PheWAS_ID,"_age"),NA)) %>%
      tidyr::drop_na(.data$age_col_complete) %>%
      dplyr::select(PheWAS_ID=.data$age_col_complete)

    all_phenos <- nearly_all_phenos %>%
      dplyr::bind_rows(age_phenos)

  } else if (!is.null(phenotype_exclusion_file)){
    if(!file.exists(phenotype_exclusion_file)){
      rlang::abort(paste0("'phenotype_exclusion_file' must be a file"))
    }
    phenotype_exclusion_ID <- data.table::fread(phenotype_exclusion_file, header = F) %>%
      dplyr::pull(1)
    nearly_all_phenos <- all_possible_phenotypes %>%
      dplyr::filter(!.data$PheWAS_ID %in% phenotype_exclusion_ID)

    age_phenos <- nearly_all_phenos %>%
      dplyr::mutate(ID2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female|_age_of_onset"),NA,.data$PheWAS_ID)) %>%
      tidyr::drop_na(.data$ID2) %>%
      dplyr::mutate(age_col_complete = ifelse(.data$age_col=="named",paste0(.data$PheWAS_ID,"_age"),NA)) %>%
      tidyr::drop_na(.data$age_col_complete) %>%
      dplyr::select(PheWAS_ID=.data$age_col_complete)

    all_phenos <- nearly_all_phenos %>%
      dplyr::bind_rows(age_phenos)
  } else {
    nearly_all_phenos <- all_possible_phenotypes

    age_phenos <- nearly_all_phenos %>%
      dplyr::mutate(ID2=ifelse(stringr::str_detect(.data$PheWAS_ID,"_male|_female|_age_of_onset"),NA,.data$PheWAS_ID)) %>%
      tidyr::drop_na(.data$ID2) %>%
      dplyr::mutate(age_col_complete = ifelse(.data$age_col=="named",paste0(.data$PheWAS_ID,"_age"),NA)) %>%
      tidyr::drop_na(.data$age_col_complete) %>%
      dplyr::select(PheWAS_ID=.data$age_col_complete)

    all_phenos <- nearly_all_phenos %>%
      dplyr::bind_rows(age_phenos)
  }

  if(!is.null(use_existing_results_file)){
    if(!file.exists(use_existing_results_file)){
      rlang::abort(paste0("'use_existing_results_file' must be a file"))
}
      old_results <- readRDS(use_existing_results_file)

      all_names <- mapply(inner_names,names(old_results),MoreArgs = list(y=old_results)) %>%
        purrr::reduce(c) %>%
        unique(.)

      all_phenos <- all_phenos %>%
        dplyr::filter(!.data$PheWAS_ID %in% all_names)

      appending_results <- TRUE
  } else {
    appending_results <- FALSE
    old_results <- NULL
    }

  if(!is.null(GRS_input)) {
      if(!file.exists(GRS_input)){
        rlang::abort(paste0("'GRS_input' must be a file"))
      }
    GRS_input_df <- data.table::fread(GRS_input)
    GRS_input_df_colname <- c("trait","group","phenotype_data","genetic_data","column_name")

      if(!all(tibble::has_name(GRS_input_df,GRS_input_df_colname))){
        warning(paste0("'GRS_input' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(GRS_input_df_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(GRS_input_df), collapse=","),
                paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(GRS_input_df),GRS_input_df_colname), collapse=","))
      }

    analysis_save_name <- NULL

    GRS_map_edit <- GRS_input_df %>%
      dplyr::mutate(save_location=paste0(analysis_folder,"/",.data$trait,"/",.data$trait,"_","results_list.RDS")) %>%
      dplyr::group_split(.data$trait)

    # create folders to host results
    mapply(making_GRS_folders,
           GRS_map_edit,
           MoreArgs = list(analysis_folder=analysis_folder))

    if(!is.na(N_cores)) {
      parallel::mcmapply(per_group_per_trait_GRS,GRS_map_edit,
                         MoreArgs = list(GRS_input=GRS_input,
                                         non_GRS_data=non_GRS_data,
                                         covariates=covariates,
                                         all_phenos=all_phenos,
                                         PheWAS_manifest=PheWAS_manifest,
                                         analysis_save_name=analysis_save_name,
                                         quant_min_cases=quant_min_cases,
                                         binary_min_cases=binary_min_cases,
                                         covariate_col_names=covariate_col_names,
                                         old_results=old_results,
                                         no_covariate=no_covariate),
                         mc.cores = N_cores)
    } else {
      mapply(per_group_per_trait_GRS,GRS_map_edit,
             MoreArgs = list(GRS_input=GRS_input,
                             non_GRS_data=non_GRS_data,
                             covariates=covariates,
                             all_phenos=all_phenos,
                             PheWAS_manifest=PheWAS_manifest,
                             analysis_save_name=analysis_save_name,
                             quant_min_cases=quant_min_cases,
                             binary_min_cases=binary_min_cases,
                             covariate_col_names=covariate_col_names,
                             old_results=old_results,
                             no_covariate=no_covariate))
    }

  } else if(!is.null(non_GRS_data)){

    genetic_data <- data.table::fread(non_GRS_data) %>%
      dplyr::rename(eid=1)

    tested_variables <- genetic_data %>%
      dplyr::select(-.data$eid) %>%
      colnames(.)

    analysis_save_name <- analysis_name

    # groups
    if(is.null(group_name_overide)){
      phenotype_group_name <-gsub("_.*", "", basename(unlist(strsplit(phenotype_files,","))))
    } else {
      phenotype_group_name <- unlist(strsplit(group_name_overide,","))
    }
    phenotype_files <- unlist(strsplit(phenotype_files,","))

    # making folders
    lapply(phenotype_group_name, function(x) dir.create(paste0(analysis_folder,"/",x),recursive = T))

    if(!is.na(N_cores)) {
      all_results <- parallel::mcmapply(per_group_none_GRS,
                                        phenotype_group_name,
                                        phenotype_files,
                                        MoreArgs = list(tested_variables=tested_variables,
                                                        GRS_input=GRS_input,
                                                        non_GRS_data=non_GRS_data,
                                                        covariates=covariates,
                                                        all_phenos=all_phenos,
                                                        PheWAS_manifest=PheWAS_manifest,
                                                        analysis_save_name=analysis_save_name,
                                                        quant_min_cases=quant_min_cases,
                                                        binary_min_cases=binary_min_cases,
                                                        covariate_col_names=covariate_col_names,
                                                        no_covariate=no_covariate),
                                        mc.cores = N_cores,
                                        SIMPLIFY = F,
                                        USE.NAMES = T)

    } else {
      all_results <-  mapply(per_group_none_GRS,
                             phenotype_group_name,
                             phenotype_files,
                             MoreArgs = list(tested_variables=tested_variables,
                                             GRS_input=GRS_input,
                                             non_GRS_data=non_GRS_data,
                                             covariates=covariates,
                                             all_phenos=all_phenos,
                                             PheWAS_manifest=PheWAS_manifest,
                                             analysis_save_name=analysis_save_name,
                                             quant_min_cases=quant_min_cases,
                                             binary_min_cases=binary_min_cases,
                                             covariate_col_names=covariate_col_names,
                                             no_covariate=no_covariate),
                             SIMPLIFY = F,
                             USE.NAMES = T)
    }


    if(is.null(old_results)){
      saveRDS(all_results,paste0(analysis_folder,"/",analysis_name,"_all_group_results_list.RDS"))
    } else {

      to_change <- old_results[names(all_results)]

      updated_results <- mapply(adding_to_results,names(to_change),d=tested_variables,e=phenotype_group_name,MoreArgs = list(b=old_results,c=all_results),SIMPLIFY = F)

      unchanged <- old_results[names(old_results)[!names(old_results) %in% names(updated_results)]]
      if(length(unchanged)>0){
        final <- append(unchanged,updated_results)[names(old_results)]
      } else {
        final <- updated_results
      }
      saveRDS(final,paste0(analysis_folder,"/",analysis_name,"_all_group_results_list.RDS"))
    }
  }
}
