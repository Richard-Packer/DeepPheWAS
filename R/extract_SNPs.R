#' Extracts SNPs from bgen or plink binary files into single file that is used for association.
#'
#' @param x Chromosome to search
#' @param analysis_folder Location of analysis
#' @param snp_guide edited version of snp_list
#' @param plink_exe command to activate plink2 in system
#' @param bgenix_exe command to activate bgenix in system
#' @param plink_type bed or pgen file if using plink files
#' @param snp_list_file location of snp_list in temp_plink file
#' @param plink_input logical value if using plink as store of genetic data
#' @param bgen_input logical value if using bgen as store of genetic data
#' @param ref_bgen option required for extracting data from begn files.
#' @return A pgen file of the selected SNPs
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
snp_selection <- function(x,analysis_folder,snp_guide,plink_exe,bgenix_exe,plink_type,snp_list_file,plink_input,bgen_input,ref_bgen) {

  output_name <- paste0(analysis_folder,"/temp_plink/",paste0(x,"_temp"))

  genetic_file_per_chromo <- snp_guide %>%
    dplyr::filter(.data$chromosome==x)

  genetic_file <- unique(genetic_file_per_chromo$genetic_file_location)

  if (bgen_input) {
    if(!dir.exists(paste0(analysis_folder,"/temp_bgen"))) {
      dir.create(paste0(analysis_folder,"/temp_bgen"),recursive = T)
    }
    temp_bgen <- paste0(analysis_folder,"/temp_bgen/",paste0(x,"_temp"))
    bgen_file <- paste0(genetic_file,".bgen")
    ref_bgen <- paste0(ref_bgen)
    sample_file <-  unique(genetic_file_per_chromo$psam_fam_sample_file_location)

    system(paste0(bgenix_exe," -g ",bgen_file," -incl-rsids ",snp_list_file," > ",temp_bgen,"\n",
                  plink_exe," --bgen ",temp_bgen," ",ref_bgen," --sample ",sample_file," --make-pgen --out ",output_name))

  } else if (plink_input) {

    if(plink_type=="bed") {
      bed_file <- paste0(genetic_file,".bed")
      bim_file <- paste0(genetic_file,".bim")
      fam_file <- unique(genetic_file_per_chromo$psam_fam_sample_file_location)

      system(paste0(plink_exe," --bed ",bed_file," --bim ",bim_file," --fam ",fam_file," --extract ",snp_list_file," --make-pgen --no-pheno --out ",output_name))

    } else if(plink_type=="pgen") {
      pgen_file <- paste0(genetic_file,".pgen")
      pvar_file <- paste0(genetic_file,".pvar")
      psam_file <- unique(genetic_file_per_chromo$psam_fam_sample_file_location)

      system(paste0(plink_exe," --pgen ",pgen_file," --pvar ",pvar_file," --psam ",psam_file," --extract ",snp_list_file," --make-pgen --multiallelics-already-joined --no-psam-pheno --out ",output_name))

    }

  }
}
#' Recodes SNP list input
#'
#' @param x SNP list input
#' @return Recoded SNPs creates unique value unlike rsid which has duplicates.
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
recoding_variants <- function(x) {

  new_pvar <- data.table::fread(x) %>%
    dplyr::mutate(REF_number=match(stringr::str_sub(.data$REF,1,1),LETTERS),
                  ALT_number=match(stringr::str_sub(.data$ALT,1,1),LETTERS),
                  multi=dplyr::case_when(nchar(.data$REF) > 1 & nchar(.data$ALT)== 1 ~ 1,
                                         nchar(.data$ALT) > 1 & nchar(.data$REF)== 1 ~ 2,
                                         nchar(.data$REF) > 1 & nchar(.data$ALT) <1 ~ 3),
                  allele_order=ifelse(.data$ALT_number>.data$REF_number,paste0(.data$REF,"_",.data$ALT),
                                      ifelse(.data$REF_number>.data$ALT_number,paste0(.data$ALT,"_",.data$REF),
                                             ifelse(.data$REF_number==.data$ALT_number,
                                                    ifelse(.data$multi==1,paste0(.data$ALT,"_",.data$REF),
                                                           ifelse(.data$multi==2,paste0(.data$REF,"_",.data$ALT),
                                                                  ifelse(nchar(.data$REF)>nchar(.data$ALT),paste0(.data$ALT,"_",.data$REF),paste0(.data$REF,"_",.data$ALT)))),NA))),
                  ID=paste0(.data$ID,"_",.data$allele_order)) %>%
    dplyr::select(.data$`#CHROM`,.data$POS,.data$ID,.data$REF,.data$ALT)

  data.table::fwrite(new_pvar,x,sep = "\t", na = NA,quote = F)

}

#' Extracts SNPs from bgen or plink binary files into single file that is used for association.
#'
#' @param genetic_file_guide Full path of the completed genetic_file_guide_template.csv.
#' @param SNP_list Full path of the file describing which SNPs are to be extracted from the genetic data files ahead of association testing.
#' @param analysis_folder Full path of the folder that will contain the data for the SNPs given by SNP_list. A temporary folder, named "temp_plink" by default, will be created within the folder specified by this argument as a place to hold temporary files needed for the SNP extraction process. Unless otherwise requested by specifying the no_delete_temp argument, the temporary folder and its contents will be deleted following successful SNP extraction.
#' @param bgen_input Specify that the genetic data files are in .bgen format.
#' @param plink_input Specify that the genetic data files are in PLINK format (.bed or .pgen).
#' @param plink_exe Full path to the PLINK2 executable.
#' @param plink_type Specify whether the PLINK-formatted genetic data are in .bed or .pgen format.
#' @param ref_bgen One of the following three values that specifies which allele is to be used as the reference: ref-first (first allele is the reference, default), ref-last (last allele is the reference), red-unknown (last allele is provisionally treated as the reference).
#' @param bgenix_exe Full path to the bgenix executable.
#' @param variant_save_name Name of the output genetic data files.
#' @param no_delete_temp Specify whether the temporary folder containing intermediate files should be retained (TRUE) or deleted (FALSE). Default is TRUE.
#' @return A pgen file of the selected SNPs
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

snp_extractor <- function(genetic_file_guide,
                          SNP_list,
                          analysis_folder,
                          bgen_input,
                          plink_input,
                          plink_exe,
                          plink_type,
                          ref_bgen,
                          bgenix_exe,
                          variant_save_name,
                          no_delete_temp) {

  # reading in the SNPs to find
  if(!file.exists(SNP_list)){
    rlang::abort(paste0("'SNP_list' must be a file"))
  }
  snp_search <- data.table::fread(SNP_list)
  snp_search_colname <- c("chromosome","rsid","group_name","coded_allele","non_coded_allele","graph_save_name")

  if(!all(tibble::has_name(snp_search,snp_search_colname))){
    warning(paste0("'SNP_list' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(snp_search_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(snp_search), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(snp_search),snp_search_colname), collapse=","))

  }

  snp_search <- snp_search %>%
    dplyr::mutate(chromosome=as.character(.data$chromosome))
  if(!file.exists(genetic_file_guide)){
    rlang::abort(paste0("'genetic_file_guide' must be a file"))
  }
  genetic_file_guide <- data.table::fread(genetic_file_guide)
  genetic_file_guide_colname <- c("chromosome","genetic_file_location","psam_fam_sample_file_location")

  if(!all(tibble::has_name(genetic_file_guide,genetic_file_guide_colname))){
    warning(paste0("'genetic_file_guide' does not have the correct colnames and may not produce the correct output, expected colnames are:
                                 '"),paste(genetic_file_guide_colname, collapse=","),paste0("'
                                 not:
                                 "),paste(colnames(genetic_file_guide), collapse=","),
            paste0("
                                 differences between inputed file and expected are:
                                 "),paste(setdiff_all(names(genetic_file_guide),genetic_file_guide_colname), collapse=","))

  }


  genetic_file_guide <- genetic_file_guide %>%
    dplyr::mutate(chromosome=as.character(.data$chromosome))

  # PLINK2 + begnix executable commands
  plink_exe <- plink_exe
  bgenix_exe <- bgenix_exe

  # folder to store the temporary Plink files
  if(!dir.exists(paste0(analysis_folder,"/temp_plink"))){
    dir.create(paste0(analysis_folder,"/temp_plink"),recursive = T)
  }

  # selects SNPs from chromosomes
  # acceptable chromosome coding
  chrom_accept <- genetic_file_guide %>%
    dplyr::select(.data$chromosome) %>%
    dplyr::pull()

  # check for acceptable chromosome input
  snp_filter <- snp_search %>%
    dplyr::filter(.data$chromosome %in% chrom_accept)

  # creates a file of SNPs in a format Plink can use to extract SNPs always goes in the same place in name analysis folder normally
  snp_list <- snp_filter %>%
    dplyr::select(.data$rsid)
  data.table::fwrite(snp_list,paste0(analysis_folder,"/temp_plink/",paste0("snp_list.txt")),col.names = F)

  # location of said file is required
  snp_list_file <- paste0(analysis_folder,"/temp_plink/",paste0("snp_list.txt"))

  snp_guide <- snp_filter %>%
    dplyr::left_join(genetic_file_guide, by="chromosome")

  chr_to_search <- unique(snp_guide$chromosome)

  # extract snps per chromosome
  mapply(snp_selection,
         chr_to_search,
         MoreArgs = list(analysis_folder=analysis_folder,
                         snp_guide=snp_guide,
                         plink_exe=plink_exe,
                         bgenix_exe=bgenix_exe,
                         plink_type=plink_type,
                         snp_list_file=snp_list_file,
                         plink_input=plink_input,
                         bgen_input=bgen_input,
                         ref_bgen=ref_bgen),
         SIMPLIFY = F, USE.NAMES = T)

  # pre_merge_psam edit to create unique ID for SNP input used in later analysis
  pvar_names <- list.files(paste0(analysis_folder,"/temp_plink"),pattern = "temp.pvar")
  pvar_file_loc <- lapply(pvar_names, function(x) paste0(analysis_folder,"/","temp_plink","/",x))

  lapply(pvar_file_loc,recoding_variants)

  # list all pgen files
  plink_merge_list <- list.files(paste0(analysis_folder,"/temp_plink"),pattern = "_temp.pgen")

  if(length(plink_merge_list)==1) {

    plink_name <- stringr::str_remove(unlist(plink_merge_list),".pgen")
    plink_location <- paste0(analysis_folder,"/temp_plink/",plink_name)
    plink_output_location <- paste0(analysis_folder,"/",variant_save_name)
    system(paste0(plink_exe," --make-pgen --pfile ",plink_location," --out ",plink_output_location))

  } else {
    #convert to df for saving
    plink_merge_edit <- data.frame(paste0(paste0(analysis_folder,"/temp_plink/",(stringr::str_remove(plink_merge_list,".pgen")))))
    data.table::fwrite(plink_merge_edit,paste0(analysis_folder,"/temp_plink/merge_list"),col.names = F)

    merge_file_location <- paste0(analysis_folder,"/temp_plink/merge_list")
    merge_output_location <- paste0(analysis_folder,"/",variant_save_name)

    #perform merge
    system(paste0(plink_exe," --pmerge-list ",merge_file_location," pfile --out ",merge_output_location))
  }

  if(no_delete_temp) {

  } else {
    unlink(paste0(analysis_folder,"/temp_plink"),recursive = T)
    if(dir.exists(paste0(analysis_folder,"/temp_bgen"))) {
      unlink(paste0(analysis_folder,"/temp_bgen"),recursive = T)
    }
  }
}
