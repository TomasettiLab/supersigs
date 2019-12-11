# TransformData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 10, 2019

# library(dplyr)
# library(here)
# source(here("code", "GenerateMinSigmaAlgebra.R"))

#' Transform data frame of mutations
#' 
#' Transform data frame of mutations by projecting counts for candidate features
#' 
#' @param dt a data frame of mutations
#' @param features_context_0 a vector of candidate features for the unexposed group
#' returned by \code{ContextMatters}
#' @param features_context_1 a vector of candidate features for the exposed group
#' returned by \code{ContextMatters}
#' @param new_partition a partition of features from \code{GenerateMinSigmaAlgebra}
#' 
#' @import dplyr
#' 
#' @return \code{TransformData} returns a transformed data frame of mutations, with
#' columns corresponding to the candidate features with projected counts
#' and other necessary columns (\code{IndVar}, \code{AGE}, \code{TOTAL_MUTATIONS},
#' and \code{DIVISON})
#' 
TransformData <- function(dt, 
                          features_context_0, 
                          features_context_1,
                          new_partition, 
                          factor){
  # Add columns for every type of mutation
  dt <- dt %>% 
    mutate_(.dots = muts_formula) %>%
    mutate(tracking_ind = 1:nrow(dt))
  
  # Store variables
  features_selected <- names(new_partition)
  
  # Count features by projections or directly summing
  if(factor != "AGE"){
    features_context_ls <- list(features_context_0, features_context_1)
    dt_new_ls <- vector("list", 2)
    
    # Project features separately for IndVar = 0 and IndVar = 1
    for(j in 1:2){
      features_context <- features_context_ls[[j]]
      # Generate partition within each IndVar group first
      features_context_partition <- GenerateMinSigmaAlgebra(input_ls = list("V1" = features_context,
                                                                            "V1" = features_context))$new_partition
      dt_new_ls[[j]] <- dt %>% 
        filter(IndVar == (j-1))
      
      # Loop over every feature
      for(i in 1:length(new_partition)){
        feature_name <- names(new_partition[i])
        feature_mutations <- new_partition[[i]]
        
        # Case 1: Feature is Remaining, then use original count
        if(feature_name == "Remaining"){
          dt_new_ls[[j]][feature_name] <- dt_new_ls[[j]][feature_mutations] %>% rowSums()
          next
        }
        
        # Case 2: Find the feature in features_context_partition 
        # that contains feature_mutations and project from it.
        # If exactly equal, this is equivalent to using the original count
        parent <- which(sapply(features_context_partition, 
                               function(f) length(setdiff(feature_mutations, f)) == 0)) %>% names()
        
        parent_prop <- sum(sapply(features_context_partition[[parent]], function(x) prop_muts_all[[x]][["TOTAL_MUTATIONS"]]))
        child_prop <- sum(sapply(feature_mutations, function(x) prop_muts_all[[x]][["TOTAL_MUTATIONS"]]))
        prop_of_parent <- child_prop/parent_prop
        dt_new_ls[[j]][feature_name] <- prop_of_parent*(dt_new_ls[[j]] %>% select(features_context_partition[[parent]]) %>% rowSums())
      }
    }
    dt_new <- bind_rows(dt_new_ls[[1]], dt_new_ls[[2]])
  } else {
    new_partition_formula <- sapply(new_partition, function(x) 
      paste0("(", paste("`", x, "`", sep = "", collapse = "+"), ")"))
    
    dt_new <- dt %>% 
      mutate_(.dots = new_partition_formula)
  }
  
  dt_new <- dt_new %>%
    arrange(tracking_ind) %>%
    select(c(features_selected, "IndVar", "AGE", "TOTAL_MUTATIONS", matches("DIVISION"))) 
  
  return(dt_new)
}

# Load data dependencies
# prop_muts_all <- readRDS(here("data", "prop_muts_all.rds"))
# muts_formula <- readRDS(here("data", "muts_formula.rds"))