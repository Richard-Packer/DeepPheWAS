% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plink_association_testing.R
\name{association_tests}
\alias{association_tests}
\title{Runs Plink2 association testing for the extracted variants}
\usage{
association_tests(
  y,
  x,
  split_groups,
  analysis_folder,
  psam,
  N_binary_split,
  N_quant_split,
  binary_ID,
  quant_ID,
  check_existing_results,
  save_plink_tables,
  analysis_name,
  plink_exe,
  association_variants,
  model,
  covariates,
  split_analysis
)
}
\arguments{
\item{y}{Group name}

\item{x}{Phenotype data frame}

\item{split_groups}{Group names(s) that are to be split}

\item{analysis_folder}{folder hosting the outputs of the function}

\item{psam}{the psam file of the extracted SNPs}

\item{N_binary_split}{numerical value for how many binary phenotypes to include in each split.}

\item{N_quant_split}{numerical value for how many quantitative phenotypes to include in each split.}

\item{binary_ID}{Binary phenotype PheWAS_IDs}

\item{quant_ID}{Quantitative phenotype PheWAS_IDs}

\item{check_existing_results}{logical value used whether to look for existing Plink results in analysis folder structure and to not repeat those associations.}

\item{save_plink_tables}{all available phenotypes generated by the previous functions.}

\item{analysis_name}{name of the analysis for saving.}

\item{plink_exe}{required input in system to run Plink2}

\item{association_variants}{extracted SNPs}

\item{model}{genetic model used in Plink}

\item{covariates}{covariates file for association analysis}

\item{split_analysis}{logical value if split analysis is being used.}
}
\value{
Plink association results.
}
\description{
Runs Plink2 association testing for the extracted variants
}
\keyword{internal}
