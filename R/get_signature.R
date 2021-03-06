# get_signature.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo, Yifan Zhang
# Date last modified: Mar 12, 2019
#
# (Export) Function to calculate signature

#' Function to obtain a SuperSig
#' 
#' Generate a tissue-specific SuperSig for a given dataset of mutations and
#' exposure factor. Returns the SuperSig and a classification model trained 
#' with the SuperSig.
#' 
#' @param data a data frame of mutations containing columns for
#' \code{sample_id}, \code{age}, \code{IndVar}, and the 96 trinucleotide
#' mutations (see vignette for details)
#' @param factor the factor/exposure (e.g. "age", "smoking"). If the 
#' factor = "age", the SuperSig is computed using
#' counts. Otherwise, rates (counts/age) are used.
#' @param wgs logical value indicating whether sequencing data is 
#' whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE}) 
#' 
#' @import dplyr
#' 
#' @export
#' 
#' @return \code{get_signature} returns an object of class \code{SuperSig}
#' 
#' @examples
#' 
#' head(example_dt) # use example data from package
#' input_dt <- make_matrix(example_dt) # convert to correct format
#' input_dt$IndVar <- c(1, 1, 1, 0, 0) # add IndVar column
#' get_signature(data = input_dt, factor = "Age") # get SuperSig
#' 
get_signature <- function(data, factor, wgs = FALSE){
    # Capitalize factor string
    factor <- toupper(factor)
    
    # Check column names of data
    if(is.na(match("SAMPLE_ID",toupper(colnames(data))))){
        stop('Input data frame missing sample_id column.')
    } else {
        colnames(data)[match("SAMPLE_ID",toupper(colnames(data)))]<-"sample_id"
    }
    if(is.na(match("AGE",toupper(colnames(data))))){
        stop('Input data frame missing AGE column.')
    } else {
        colnames(data)[match("AGE",toupper(colnames(data)))] <- "AGE"
    }
    if(is.na(match("INDVAR",toupper(colnames(data))))){
        stop('Input data frame missing IndVar column.')
    } else {
        colnames(data)[match("INDVAR",toupper(colnames(data)))] <- "IndVar"
    }

    if(nrow(data) < 5) stop('More than 5 samples required for get_signature.')
    
    # Check all trinucleotides are present
    trinucleotideBases <- unique(transform_muts_vec)
    if(all(trinucleotideBases %in% toupper(colnames(data)))){
        data$TOTAL_MUTATIONS <- rowSums(data[,trinucleotideBases])
    } else {
        stop('Input data frame missing one or more trinucleotide mutations.')
    }
    
    # Get features
    features_out <- suppressWarnings(feature_selection(data, factor, wgs))
    
    # Get apparent AUC and model
    classification_out <- supersig_classifier(dt = features_out$dt_new,
                                              factor = factor,
                                              keep_classifier = TRUE,
                                              features_selected = 
                                                features_out$features_selected,
                                              select_n = features_out$select_n)
    
    # Create and return S4 object
    feature_names <- names(classification_out$signature$mean_diffs)
    out <- SuperSig(Signature = classification_out$signature$mean_diffs,
                    Features = features_out$new_partition[feature_names],
                    AUC = classification_out$auc,
                    Model = classification_out$classifier)
    return(out)
}