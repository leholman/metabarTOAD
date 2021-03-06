% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DenoiseUNOISE3.R
\name{DenoiseUNOISE3}
\alias{DenoiseUNOISE3}
\title{DenoiseUNOISE3}
\usage{
DenoiseUNOISE3(UsePrimerFile = FALSE, usearchdest = "usearch")
}
\arguments{
\item{UsePrimerFile}{TRUE or FALSE, is a csv file provided containing primer information?}

\item{usearchdest}{Specify the location of the usearch executable if PATH not configured appropriately.}
}
\value{
None
}
\description{
This function denoises amplicons using the UNOISE3 algorithm from the USEARCH softwae by Robert Edgar.
It depends on the folder structure created by the function Folders().
}
