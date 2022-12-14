% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_filtering.R
\name{tab_data_clean}
\alias{tab_data_clean}
\title{Searches and extracts named columns from data file}
\usage{
tab_data_clean(x, y, z)
}
\arguments{
\item{x}{full file path to data file}

\item{y}{vector of column names to extract}

\item{z}{dataframe of exclusions single column no header}
}
\value{
A dataframe containing only columns of interest
}
\description{
Searches and extracts named columns from data file
}
\keyword{internal}
