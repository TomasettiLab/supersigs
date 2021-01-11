# get_signature.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo, Yifan Zhang
# Date last modified: Mar 12, 2019
#
# (Export) Function to calculate signature

# library(dplyr)

#' Function to obtain a SuperSig
#' 
#' Generate a tissue-specific SuperSig for a given dataset of mutations and
#' exposure factor. Returns the SuperSig and a classification model trained with
#' the SuperSig.
#' 
#' @param dt a data frame of mutations (see vignette for details)
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' @param wgs logical value indicating whether sequencing data is 
#' whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE}) 
#' 
#' @import dplyr
#' 
#' @export
#' 
#' @return \code{get_signature} returns an object of class "SuperSig"
#' 
#' @examples
#' 
#' # print(example_dt) # use example data from package
#' # input_dt = make_matrix(example_dt) # convert to correct format
#' # get_signature(dt = input_dt, factor = "age") # get SuperSig
#' 
get_signature <- function(dt, factor, wgs = F){
  # Capitalize factor string
  factor = toupper(factor)
  
  # Check column names of dt
  if(is.na(match("SAMPLE_ID",toupper(colnames(dt))))){
    stop('Input data frame missing sample_id column.')
  } else {
    colnames(dt)[match("SAMPLE_ID",toupper(colnames(dt)))] = "sample_id"
  }
  if(is.na(match("AGE",toupper(colnames(dt))))){
    stop('Input data frame missing AGE column.')
  } else {
    colnames(dt)[match("AGE",toupper(colnames(dt)))] = "AGE"
  }
  if(is.na(match("INDVAR",toupper(colnames(dt))))){
    stop('Input data frame missing IndVar column.')
  } else {
    colnames(dt)[match("INDVAR",toupper(colnames(dt)))] = "IndVar"
  }
  
  # Check number of samples (> 5 required)
  if(nrow(dt) < 5) {
    stop("More than 5 samples are required to run get_signature.")
  }
  
  # Check if counts of 96 trinucleotide bases are present, and compute total mutations
  trinucleotideBases <- unique(transform_muts_vec)
  if(all(trinucleotideBases %in% toupper(colnames(dt)))){
    dt$TOTAL_MUTATIONS <- rowSums(dt[,trinucleotideBases])
  } else {
    stop('Input data frame missing one or more trinucleotide mutations.')
  }
  
  # Get features
  features_out = suppressWarnings(feature_selection(dt, factor, wgs))
  
  # Get apparent AUC and model
  classification_out = supersig_classifier(dt = features_out$dt_new,
                                          factor = factor,
                                          classifier = c("Logit"),
                                          keep_classifier = T,
                                          features_selected = features_out$features_selected,
                                          select_n = features_out$select_n)
  
  # Create S4 object output
  out = SuperSig(Signature = classification_out$signature$mean_diffs,
                  Features = features_out$new_partition[names(classification_out$signature$mean_diffs)],
                  AUC = classification_out$auc,
                  Model = classification_out$classifier)
  
  return(out)
}