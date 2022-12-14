% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primary_care_quant_phenotypes.R
\name{quant_primary_care}
\alias{quant_primary_care}
\title{Makes primary care quantitative phenotypes guided by PheWAS manifest}
\usage{
quant_primary_care(a, b, c, d, e, prim_care, DOB, code_lists)
}
\arguments{
\item{a}{PheWAS_ID}

\item{b}{primary care code list}

\item{c}{existence of limits or not}

\item{d}{lower_limit quant value}

\item{e}{upper_limit quant value}

\item{prim_care}{primary care data}

\item{DOB}{DOB information}

\item{code_lists}{Full file path of the folder containing alternative primary care code lists.}
}
\value{
A list of primary care quantitative phenotypes
}
\description{
Makes primary care quantitative phenotypes guided by PheWAS manifest
}
\keyword{internal}
