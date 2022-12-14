% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenotype_preparation.R
\name{converted_for_association}
\alias{converted_for_association}
\title{Converts phenotypes to form suitable for association analysis retains list format}
\usage{
converted_for_association(x, y, PheWAS_manifest, PheWAS_ID_list)
}
\arguments{
\item{x}{Phenotypes as list object}

\item{y}{PheWAS_ID of phenotype}

\item{PheWAS_manifest}{the PheWAS manifest}

\item{PheWAS_ID_list}{all PheWAS ID's in analysis}
}
\value{
Phenotypes as list edited to for combination into dataframe for association analysis
}
\description{
Converts phenotypes to form suitable for association analysis retains list format
}
\keyword{internal}
