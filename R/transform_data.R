# transform_data.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Feb 15, 2021

# library(dplyr)
# library(here)

#' Transform data frame of mutations
#' 
#' Transform data frame of mutations by projecting counts for candidate features
#' 
#' @param dt a data frame of mutations
#' @param new_partition a partition of features from 
#' \code{generate_min_sigma_algebra}
#' 
#' @import dplyr
#' @importFrom rlang .data
#' 
#' @return \code{transform_data} returns a transformed data frame of mutations, 
#' with columns corresponding to the candidate features with projected counts
#' and other necessary columns (\code{IndVar}, \code{AGE}, 
#' \code{TOTAL_MUTATIONS}, and \code{DIVISON})
#' 
#' @noRd
#' 
transform_data <- function(dt, 
                           new_partition){
    # Add columns for every type of mutation
    dt <- dt %>% 
        mutate_(.dots = muts_formula) %>%
        mutate(tracking_ind = seq_len(nrow(dt)))
    
    # Store variables
    features_selected <- names(new_partition)
    
    # Count features by directly summing
    new_partition_formula <- vapply(new_partition, function(x) 
        paste0("(", paste("`", x, "`", sep = "", collapse = "+"), ")"),
        FUN.VALUE = character(1))
    
    dt_new <- dt %>% 
        mutate_(.dots = new_partition_formula)
    
    dt_new <- dt_new %>%
        arrange(.data$tracking_ind) %>%
        select(c(features_selected, "IndVar", "AGE", "TOTAL_MUTATIONS", 
                         matches("DIVISION"))) 
    
    return(dt_new)
}

# Load data dependencies
# prop_muts_all <- readRDS(here("data", "prop_muts_all.rds"))
# muts_formula <- readRDS(here("data", "muts_formula.rds"))
