# TransformData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 21, 2020

# library(dplyr)
# library(here)
# source(here("code", "GenerateMinSigmaAlgebra.R"))

#' Transform data frame of mutations
#' 
#' Transform data frame of mutations by projecting counts for candidate features
#' 
#' @param dt a data frame of mutations
#' @param new_partition a partition of features from \code{GenerateMinSigmaAlgebra}
#' 
#' @import dplyr
#' 
#' @return \code{TransformData} returns a transformed data frame of mutations, with
#' columns corresponding to the candidate features with projected counts
#' and other necessary columns (\code{IndVar}, \code{AGE}, \code{TOTAL_MUTATIONS},
#' and \code{DIVISON})
#' 
#' @noRd
#' 
TransformData <- function(dt, 
                          new_partition){
  # Add columns for every type of mutation
  dt <- dt %>% 
    mutate_(.dots = muts_formula) %>%
    mutate(tracking_ind = 1:nrow(dt))
  
  # Store variables
  features_selected <- names(new_partition)
  
  # Count features by directly summing
  new_partition_formula <- sapply(new_partition, function(x) 
    paste0("(", paste("`", x, "`", sep = "", collapse = "+"), ")"))
  
  dt_new <- dt %>% 
    mutate_(.dots = new_partition_formula)
  
  dt_new <- dt_new %>%
    arrange(tracking_ind) %>%
    select(c(features_selected, "IndVar", "AGE", "TOTAL_MUTATIONS", matches("DIVISION"))) 
  
  return(dt_new)
}

# Load data dependencies
# prop_muts_all <- readRDS(here("data", "prop_muts_all.rds"))
# muts_formula <- readRDS(here("data", "muts_formula.rds"))