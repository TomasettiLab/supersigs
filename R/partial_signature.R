# partial_signature.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Feb 25, 2021
#
# (Export) Function to remove signature contribution from data

# library(dplyr)

#' Function to remove the contribution of a SuperSig
#' 
#' Remove the contribution of a SuperSig from the data and return the data.
#' 
#' @param data a data frame of mutations (see vignette for details)
#' @param object an object of class "SuperSig"
#' 
#' @import dplyr
#' 
#' @export
#' 
#' @return \code{predict_signature} returns the original data frame with 
#' the contribution of a supervised signature removed
#' 
#' @examples
#' 
#' head(example_dt) # use example data from package
#' input_dt <- make_matrix(example_dt) # convert to correct format
#' input_dt$IndVar <- c(1, 1, 1, 0, 0) # add IndVar column
#' supersig <- get_signature(data = input_dt, factor = "age") # get SuperSig
#' partial_signature(data = input_dt, object = supersig)
#' 
partial_signature <- function(data, object){
    # Check column names of data
    if(is.na(match("SAMPLE_ID",toupper(colnames(data))))){
        stop('Input data frame missing sample_id column.')
    } else {
        colnames(data)[match("SAMPLE_ID",toupper(colnames(data)))] = "sample_id"
    }
    if(is.na(match("AGE",toupper(colnames(data))))){
        stop('Input data frame missing AGE column.')
    } else {
        colnames(data)[match("AGE",toupper(colnames(data)))] = "AGE"
    }
    if(is.na(match("INDVAR",toupper(colnames(data))))){
        stop('Input data frame missing IndVar column.')
    } else {
        colnames(data)[match("INDVAR",toupper(colnames(data)))] = "IndVar"
    }
    
    # Extract slots from object
    features <- Features(object)
    
    # Subtract contributions (rates*age)
    for(i in seq_along(features)){
        trinucleotide_features <- unname(features[[i]])
        diff <- Signature(object)[[i]]
        data <- data %>%
            mutate_at(vars(all_of(trinucleotide_features)), ~.- diff*AGE) %>%
            mutate_at(vars(all_of(trinucleotide_features)), 
                      ~ifelse(. < 0, 0, .))
    }
    
    return(data)
}