# GenerateMinSigmaAlgebra.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019
#
# Function to create the smallest non-overlapping partition that generates 
# the minimal sigma algebra containing two sets of mutation features

# library(dplyr)
# library(assertthat)
# library(here)

# source(here("code", "condense_mutations.R"))
# source(here("code", "convert_to_level3.R"))

#' Create the smallest non-overlapping partition that generates 
#' the minimal sigma algebra containing two sets of mutation features
#' 
#' @param input_ls a named list of length 2, where each element is a vector of mutation features
#' @param condense boolean toggle for whether the output new_partition should return
#' mutations in condensed form (default = FALSE)
#' 
#' @import dplyr
#' @import assertthat
#' 
#' @return output is a list of 1 element, new_partition, which is the partition of features
#' 
GenerateMinSigmaAlgebra <- function(input_ls, 
                                    condense = F){
  # Check input
  assert_that(length(input_ls) == 2, msg = "input_ls has length not equal to 2")
  assert_that(!is.null(names(input_ls)), msg = "input_ls is not named")

  # Convert to fundamental (i.e. level 3) mutations
  feat_ls <- unique(unlist(input_ls))
  feat_ls <- lapply(feat_ls, convert_to_level3)
  
  # Create matrix of indicator values
  temp <- setNames(rep(F, length(muts_level3)), muts_level3)
  M <- sapply(feat_ls, FUN = function(mutation) {out <- temp; out[mutation] <- T; out})
  
  # Separate mutations that do not appear in any of feat_ls
  any_ind <- apply(M, MARGIN = 1, FUN = any)
  other_muts <- rownames(M)[!any_ind] 
  M <- M[any_ind, ]
  
  # Create partition
  if(is.null(dim(M))){ # If there is only one feature
    new_partition <- list(names(M))
  } else {
    new_partition <- lapply(rownames(M), FUN = function(mutation){
      apply(cbind(M[, which(M[mutation, ] == T)],
                  !M[, which(M[mutation, ] == F)]),
            MARGIN = 1, FUN = all)} %>% 
        which() %>% names())
    new_partition <- new_partition[!duplicated(new_partition)]
  }
  names(new_partition) <- paste0("F", seq_along(new_partition))
  
  # Add back other mutations as "Remaining" feature
  if(length(other_muts) > 0)
    new_partition[["Remaining"]] <- other_muts 
  
  # Condense mutations
  if(condense)
    new_partition <- lapply(new_partition, condense_mutations) 
  
  # Return output
  out <- list(new_partition = new_partition)
  
  return(out)
}

# Load data dependencies
# muts_level3 = readRDS(here("data", "muts_level3.rds"))
# muts_children_level3 = readRDS(here("data", "muts_children_level3.rds"))

# Test function
# h_mix = readRDS(here("data", "h_mix.rds"))
# test_1 <- h_mix[c("A[T>G]G","C[T>G]G")]
# GenerateMinSigmaAlgebra(test_1)
# 
# test_2 <- h_mix[c("(TG)[T>C](CTG)","(AG)[T>C](CTG)")]
# GenerateMinSigmaAlgebra(test_2)
