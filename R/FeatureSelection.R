# FeatureSelection.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019
#
# Function for selecting features (formerly part of MyLDAEnvClassifier.R)

# library(dplyr)
# library(assertthat)
# library(here)
# source(here("code", "ContextMatters.R"))
# source(here("code", "GenerateMinSigmaAlgebra.R"))
# source(here("code", "TransformData.R"))
# source(here("code", "PredictiveFeatures.R"))

#' Function for feature selection
#' 
#' Perform feature selection given a dataset of mutations by calling 
#' a series of other functions to find significant and predictive features
#' 
#' @param dt a data frame of mutations
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' @param test_ind an optional vector of indices for the test data
#' @param middle_dt a data frame of the middle aged samples (\code{NULL}
#' if the factor is not age)
#' @param keep_nonpredictive logical value indicating whether to combine 
#' non-predictive features as one feature (default is \code{FALSE})
#' 
#' @import dplyr
#' @import assertthat
#' 
#' @export
#' 
#' @return \code{FeatureSelection} returns a list of several elements:
#' \itemize{
#' \item \code{features_context_0} is a vector of survival mutations
#' for the unexposed group
#' \item \code{features_context_1} is a vector of survival mutations
#' for the exposed group
#' \item \code{features_selected} is a vector of candidate features ranked by AUC
#' \item \code{select_n} is the number of top features to retain for each method
#' \item \code{dt_new} is the transformed data from \code{TransformData}
#' }
#' 
FeatureSelection <- function(dt,
                             factor,
                             test_ind = NULL, 
                             middle_dt = NULL,
                             keep_nonpredictive = F){
  if(is.null(test_ind)){
    train_ind <- 1:nrow(dt)
    test_ind <- train_ind
  } else {
    train_ind <- setdiff(1:nrow(dt), test_ind)
  }
  
  # Separate into exposed (train_1) and unexposed (train_0)
  train <- dt[train_ind, ]
  if(factor != "AGE"){ 
    train_0 <- train %>% filter(IndVar == 0)
    train_1 <- train %>% filter(IndVar == 1)
  } else {
    train_0 <- rbind(train, middle_dt) 
    train_1 <- train_0
  }
  
  # Add up counts for every mutation
  train_0 <- train_0 %>%
    transmute_(.dots = muts_formula) %>% 
    mutate(TOTAL_MUTATIONS = select(., 1:6) %>% rowSums())
  
  train_1 <- train_1 %>%
    transmute_(.dots = muts_formula) %>%
    mutate(TOTAL_MUTATIONS = select(., 1:6) %>% rowSums())
  
  # (Note: Removed Mutrelative2TotalPval section)
  
  # Test for significant features
  features_context_0 <- ContextMatters(train_0)
  assert_that(length(features_context_0) >= 1, 
              msg = "No significant features found for ind0 by ContextMatters")
  features_context_1 <- ContextMatters(train_1)
  assert_that(length(features_context_1) >= 1, 
              msg = "No significant features found for ind1 by ContextMatters")
  input_ls <- list(var0 = features_context_0, var1 = features_context_1)
  features_gmsa <- GenerateMinSigmaAlgebra(input_ls)
  new_partition <- features_gmsa$new_partition
  features_selected <- names(new_partition)
  
  # Transform data
  dt_new <- TransformData(dt,
                          features_context_0, 
                          features_context_1,
                          new_partition, 
                          factor)
  
  # Test for predictive features
  predictive_out <- PredictiveFeatures(train = dt_new[train_ind, ], 
                                       new_partition = new_partition,
                                       factor = factor)
  features_selected <- predictive_out$features_selected
  
  # Add combined non-predictive features (Rest) to features_selected
  if(keep_nonpredictive){
    dt_new <- dt_new %>%
      mutate(SumPredictive = dt_new %>% select(features_selected) %>% rowSums(),
             Rest = TOTAL_MUTATIONS - SumPredictive)
    if(length(unique(dt_new$Rest)) != 1){ # Add Rest only if Rest is not a constant vector of 0's
      features_selected <- c(features_selected, "Rest")
    }
  }
  
  # Return output
  assert_that(length(features_selected) >= 1, 
              msg = "No predictive features found")
  
  out <- list(features_context_0 = features_context_0,
              features_context_1 = features_context_1,
              features_selected = features_selected,
              select_n = predictive_out$select_n,
              dt_new = dt_new) 
  
  return(out)
}


# Load data dependencies
# muts_formula <- readRDS(here("data", "muts_formula.rds"))

# Test function
# signature_caf <- readRDS(here("super_sigs", "data", "signature_caf.rds"))
# factor <- "AGE"
# tissue <- "UCEC"
# ind <- which((signature_caf["Factor",] == factor) & (signature_caf["Tissue",] == tissue))
# dt <- signature_caf[["Data", ind]]$DataSetFiltered %>%
#   filter(TOTAL_MUTATIONS > 0)
# 
# test_ind = NULL
# middle_dt = NULL
# 
# if(factor == "AGE"){
#   middle_dt <- signature_caf[["Data", ind]]$DataSetFilteredKeepMiddle %>%
#     filter(!(SAMPLE %in% dt$SAMPLE)) %>%
#     mutate(IndVar = NA)
# }
# 
# test_out = FeatureSelection(dt = dt, middle_dt = middle_dt, factor = factor)
