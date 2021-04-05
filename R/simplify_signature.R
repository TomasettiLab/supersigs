# simplify_signature.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Mar 31, 2021
#
# (Export) Function to simplify signature representation

# Convert single feature to simplified labels (helper function)
features_to_labels <- function(feature_name, features_partition){
  mutations <- features_partition[[feature_name]]
  label_mutations <- h_mix[vapply(h_mix, function(x){setequal(x, mutations)},
                                  logical(1))]
  label <- names(label_mutations)
  if(length(label) == 1){
    return(list(label = label, 
                multiplier = 1))
  } else {
    # If no exact match in h_mix, pick labels recursively to build up feature
    label <- c()
    multiplier <- c()
    denominator <- sum(vapply(mutations, function(x){
      prop_muts_all[[x]][["TOTAL_MUTATIONS"]]},
      numeric(1)))
    while(length(mutations) > 0){
      # Pick h_mix mutation that explains as much as the feature as possible
      candidate_labels <- h_mix[vapply(h_mix, function(x){
        length(setdiff(x, mutations)) == 0}, logical(1))]
      candidate_labels_lengths <- vapply(candidate_labels, 
                                         function(x) length(x),
                                         numeric(1))
      max_candidate <- vapply(candidate_labels, function(x){
        length(x) == max(candidate_labels_lengths)}, logical(1))
      label_mutations <- candidate_labels[max_candidate][1]
      label <- c(label, names(label_mutations))
      
      # Calculate the proportion of the feature that is
      # represented by the h_mix mutation 
      # (i.e. projection based on background ratio)
      numerator <- sum(vapply(label_mutations[[1]], function(x){
        prop_muts_all[[x]][["TOTAL_MUTATIONS"]]}, numeric(1)))
      multiplier <- c(multiplier, numerator/denominator) 
      mutations <- setdiff(mutations, label_mutations[[1]])
    }
  }
  
  return(list(label = label, 
              multiplier = multiplier))
}

# Convert label to IUPAC naming
transform_to_iupac <- function(z){
  zp <- z
  for(i in seq_along(transform_iupac_grep)){
    zp <- vapply(zp, FUN = gsub,
                 pattern = transform_iupac_grep[[i]],
                 replacement = names(transform_iupac_grep)[i],
                 character(1))
  }
  return(zp[[1]])
}

#' Function to simplify signature representation 
#' into interpretable labels for visualization purposes
#' 
#' Take a signature representation from SuperSig
#' and group trinucleotides within each feature into
#' interpretable labels, with optional IUPAC labeling
#' from IUPAC_CODE_MAP in the Biostrings package
#' 
#' @param object an object of class \code{SuperSig}
#' @param iupac logical value indicating whether to use IUPAC labels
#' (iupac = \code{TRUE}) or not (iupac = \code{FALSE}) 
#' 
#' @export
#' 
#' @return \code{simplify_signature} returns a vector of
#' simplified features and their difference in mean
#' mean rates between exposed and unexposed (or
#' average rate if the factor is "age")
#' 
#' @examples
#' 
#' head(example_dt) # use example data from package
#' input_dt <- make_matrix(example_dt) # convert to correct format
#' input_dt$IndVar <- c(1, 1, 1, 0, 0) # add IndVar column
#' set.seed(1)
#' supersig <- get_signature(data = input_dt, factor = "age") # get SuperSig
#' simplify_signature(object = supersig, iupac = FALSE)
#' 
simplify_signature <- function(object, iupac){
  signature <- Signature(object)
  features_partition <- Features(object)
  signature_labeled = c()
  for(j in seq_along(signature)){
    feature_name <- names(signature[j])
    rate <- signature[[j]]
    labels_out <- features_to_labels(feature_name, features_partition)
    tmp_vec <- labels_out$multiplier*rate
    names(tmp_vec) <- labels_out$label
    signature_labeled <- c(signature_labeled, tmp_vec)
  }
  
  if(iupac){
    labels_iupac <- vapply(names(signature_labeled), 
                           function(x) transform_to_iupac(x),
                           character(1))
    names(signature_labeled) <- labels_iupac
  } 
  
  return(signature_labeled)
}
