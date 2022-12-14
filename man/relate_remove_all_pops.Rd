% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenotype_preparation.R
\name{relate_remove_all_pops}
\alias{relate_remove_all_pops}
\title{Removes related individuals from phenotypes per grouping.}
\usage{
relate_remove_all_pops(x, call_rate_kinships, PheWAS_manifest)
}
\arguments{
\item{x}{all phenotypes as list.}

\item{call_rate_kinships}{call_rate_kinship file}

\item{PheWAS_manifest}{manifest file}
}
\value{
Phenotypes filtered for relatedness.
}
\description{
Removes related individuals from phenotypes per grouping.
}
\keyword{internal}
