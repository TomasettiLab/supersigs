% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_signature.R
\name{get_signature}
\alias{get_signature}
\title{Function to obtain a SuperSig}
\usage{
get_signature(data, factor, wgs = FALSE)
}
\arguments{
\item{data}{a data frame of mutations containing columns for
\code{sample_id}, \code{age}, \code{IndVar}, and the 96 trinucleotide
mutations (see vignette for details)}

\item{factor}{the factor/exposure (e.g. "age", "smoking"). If the 
factor = "age", the SuperSig is computed using
counts. Otherwise, rates (counts/age) are used.}

\item{wgs}{logical value indicating whether sequencing data is 
whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE})}
}
\value{
\code{get_signature} returns an object of class \code{SuperSig}
}
\description{
Generate a tissue-specific SuperSig for a given dataset of mutations and
exposure factor. Returns the SuperSig and a classification model trained 
with the SuperSig.
}
\examples{

head(example_dt) # use example data from package
input_dt <- make_matrix(example_dt) # convert to correct format
input_dt$IndVar <- c(1, 1, 1, 0, 0) # add IndVar column
get_signature(data = input_dt, factor = "Age") # get SuperSig

}
