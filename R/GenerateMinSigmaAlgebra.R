# GenerateMinSigmaAlgebra.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Jul 28, 2019
#
# Function to create the smallest non-overlapping partition that generates 
# the minimal sigma algebra containing two sets of mutation features, and
# complementary helper functions convert_to_level3 and condense_mutations

library(dplyr)
library(assertthat)
library(here)

source(here("code", "condense_mutations.R"))
source(here("code", "convert_to_level3.R"))

# in$input_ls is a named list of length 2. Each element is a vector of mutation features
# in$condense toggles whether new_partition should return mutations in condensed form
# out$new_partition is the partition (note: new_partition_formula has been moved to TransformData)
# out$sup_partition is formula of input_ls based on the new partition (mainly for testing purposes)
# out$sup_partition_formula is the corresponding formula string
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
  if(length(other_muts) > 0)
    new_partition[["Remaining"]] <- other_muts # add back the other mutations
  if(condense)
    new_partition <- lapply(new_partition, condense_mutations) 
  
  # Get the partition values that make up the original input feat_ls
  sup_partition <- sapply(input_ls, FUN = function(feats){
    sapply(new_partition, FUN = function(feats_new){
      length(intersect(feats, feats_new)) > 0})}) %>%
    apply(MARGIN = 2, FUN = function(col) names(which(col)))
  
  if(!is.null(dim(sup_partition))){
    sp_names <- colnames(sup_partition)
    sup_partition <- lapply(seq_len(ncol(sup_partition)), function(i) sup_partition[,i])
    names(sup_partition) <- sp_names
  }
  
  # Return output
  out <- list(new_partition = new_partition,
              sup_partition = sup_partition,
              sup_partition_formula = sapply(sup_partition, function(x) 
                paste("`", x, "`", sep = "", collapse = "+")))
  
  return(out)
}

# Load data dependencies
muts_level3 = readRDS(here("data", "muts_level3.rds"))
muts_children_level3 = readRDS(here("data", "muts_children_level3.rds"))

# Test function
# h_mix = readRDS(here("data", "h_mix.rds"))
# test_1 <- h_mix[c("A[T>G]G","C[T>G]G")]
# GenerateMinSigmaAlgebra(test_1)
# 
# test_2 <- h_mix[c("(TG)[T>C](CTG)","(AG)[T>C](CTG)")]
# GenerateMinSigmaAlgebra(test_2)
