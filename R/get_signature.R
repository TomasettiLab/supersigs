# get_signature.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 11, 2019
#
# (Export) Function to calculate signature

# library(dplyr)
library(here)

#' Function to obtain a SuperSig
#' 
#' Generate a tissue-specific SuperSig for a given dataset of mutations and
#' exposure factor. Returns the SuperSig and a classification model trained with
#' the SuperSig.
#' 
#' @param dt a data frame of mutations
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' 
#' @import dplyr
#' 
#' @export
#' 
#' @return \code{get_signature} returns a list of several elements:
#' \itemize{
#' \item \code{Signature} is a vector of the signature, i.e. the features that
#' comprise the signature and their difference in mean rates (or counts if the factor
#' is "age")
#' \item \code{Features} is a list of features that comprise the signature and 
#' their representation in terms of the fundamental (trinucleotide) mutations
#' \item \code{AUC} is the apparent AUC (i.e. not cross-validated) obtained within 
#' the provided data using the generated SuperSig
#' \item \code{Model} is the logistic regression model trained on the data and
#' the generated SuperSig
#' }
#' 
#' @examples
#' get_signature(mutation_data, "age")
#' 
get_signature <- function(dt, factor){
  
  # Check dt argument
  if (is.na(match("SAMPLE_ID",toupper(colnames(dt))))) {
    stop('Data frame missing sample_id column.') } else {
      colnames(dt)[match("SAMPLE_ID",toupper(colnames(dt)))]="sample_id"
    }
  if (is.na(match("AGE",toupper(colnames(dt))))) {
    stop('Data frame missing AGE column.') } else {
      colnames(dt)[match("AGE",toupper(colnames(dt)))]="AGE"
    }
  if (is.na(match("INDVAR",toupper(colnames(dt))))) {
    stop('Data frame missing IndVar column.') } else {
      colnames(dt)[match("INDVAR",toupper(colnames(dt)))]="IndVar"
    }
  
  # Check if counts of 96 trinucleotide bases are present, and compute total mutations
  trinucleotideBases <- unique(transform_muts_vec)
  if (all(trinucleotideBases %in% toupper(colnames(dt)))) {
    dt$TOTAL_MUTATIONS <- rowSums(dt[,trinucleotideBases],2)
  } else {
    stop('There needs to be 96 trinucleotide mutation columns.')
  }
  
  # Get features
  features_out = suppressWarnings(FeatureSelection(dt, factor))
  
  # Get apparent AUC and model
  classification_out = SuperSigClassifier(dt = features_out$dt_new,
                                          factor = factor,
                                          classifier = c("Logit"),
                                          keep_classifier = T,
                                          features_selected = features_out$features_selected,
                                          select_n = features_out$select_n)
  
  out = list(Signature = classification_out$signature$mean_diffs,
             Features = features_out$new_partition[names(classification_out$signature$mean_diffs)],
             AUC = classification_out$auc,
             Model = classification_out$classifier)
  
  return(out)
}