# PredictiveFeatures.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019
#
# Function for testing and choosing predictive features (formerly part of MyLDAEnvClassifier.R)

# library(dplyr)
# library(caret)
# library(rsample)
# library(here)
# source(here("code", "MyAuc.R"))

#' Select predictive features
#' 
#' Rank candidate features by their median bootstrap AUC and find the top n features
#' for each classifier method, where n is chosen by cross-validation
#' 
#' @param train a transformed data of mutations from \code{TransformData}
#' @param new_partition a partition of features from \code{GenerateMinSigmaAlgebra}
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' @param n_iter an optional parameter specifying the number of 
#' cross-validation iterations (default is \code{5})
#' @param n_fold an optional parameter specifying the number of folds
#' for each cross-validation iteration (default is \code{2})
#' 
#' @import dplyr
#' @import caret
#' @import rsample
#' 
#' @export
#' 
#' @return `PredictiveFeatures` returns a list of two elements:
#' \itemize{
#' \item \code{features_selected} is a vector of candidate features ranked by AUC
#' \item \code{select_n} is the number of top features to retain for each method,
#' chosen by cross-validation
#' }
#' 
PredictiveFeatures <- function(train, 
                               new_partition,
                               factor,
                               n_iter = 5,
                               n_fold = 2){
  # Use rates for non-age factors
  if(factor != "AGE"){
    train <- train %>%
      mutate_at(vars(names(new_partition)), ~(./AGE))
  }
  
  # Initialize variables
  feature_names <- names(new_partition)
  
  # Bootstrap AUC
  boot_dt <- train %>%
    select(c(feature_names, "IndVar")) %>%
    bootstraps(1000, strata = "IndVar")
  
  auc_mat <- sapply(boot_dt$splits, function(x) 
    AllMyAuc(as_tibble(x), IndVar = "IndVar"))
  
  if(factor == "AGE"){
    auc_mat <- auc_mat  # Only count AUC > 0.5 for age
  } else {
    auc_mat <- abs(auc_mat - 0.5) + 0.5 # Count AUC > 0.5 or < 0.5 for other variables
  }
  
  # Compute median AUC
  auc_medians <- apply(auc_mat, 1, function(x) median(x, na.rm = T))
  auc_medians <- sapply(auc_medians, function(x) ifelse(is.na(x), 0, x)) # Replace missing AUCs with 0
  
  # Create mutation_dt
  mutation_dt <- tibble(mutation = feature_names,
                        auc = auc_medians)
  
  # For age, throw out candidate features with AUC medians less than 0.5 or negative differences in means
  if(factor == "AGE"){
    feature_means <- train %>%
      select(c(feature_names, "IndVar")) %>%
      group_by(IndVar) %>%
      summarize_all(mean)
    
    feature_means_diff <- feature_means %>% filter(IndVar == TRUE) - 
      feature_means %>% filter(IndVar == FALSE)
    
    remove_neg_features <- feature_means_diff %>%
      select_if(function(col) col < 0) %>%
      names()
    
    mutation_dt <- mutation_dt %>%
      filter(auc > 0.5) %>%
      filter(!(mutation %in% remove_neg_features))
  }
  
  # Inner loop of CV to find the best select_n for each classifier
  max_n = nrow(mutation_dt)
  auc_by_methods_mean = vector("list", max_n)
  methods <- c("LDA", "Logit", "RF")
  print(paste("Begin cross-validated selection over", max_n, "features..."))
  for(k in 1:max_n){
    print(paste("...testing n_features =", k))
    inner_partitions <- sapply(1:n_iter, FUN = function(x) createFolds(y = train$IndVar, k = n_fold))
    inner_aucs <- vector(length = length(inner_partitions), mode = "list")
    
    for(ij in seq_along(inner_partitions)){
      i <- (ij-1) %% n_fold + 1
      j <- floor((ij-1)/n_fold) + 1
      
      features_selected <- mutation_dt %>%
        top_n(n = k, wt = auc) %>%
        pull(mutation)
      
      inner_aucs[[ij]] = MyLDAEnvClassifier(dt = train, 
                         test_ind = inner_partitions[[i, j]],
                         factor,
                         classifier = methods, 
                         keep_classifier = F,
                         adjusted_formula = F,
                         features_selected = features_selected,
                         select_n = c("LDA" = k, "RF" = k, "Logit" = k))
    }
    
    inner_aucs_mat <- matrix(inner_aucs, nrow = n_fold)
    auc_by_methods <- sapply(methods, function(method) 
      apply(inner_aucs_mat, MARGIN = c(1, 2), function(x){x[[1]]$auc[method]}), simplify = F)
    auc_by_methods_mean[[k]] <- sapply(auc_by_methods, function(x) mean(x, na.rm = T))
  }
  
  # Best select_n for each method
  select_n = sapply(methods, function(classifier){
    tmp_max = which.max(unlist(lapply(auc_by_methods_mean, function(auc) auc[classifier])))
    if(identical(tmp_max, integer(0))) return(NA) # If classifier fails for all inner folds
    else return(tmp_max)})
  names(select_n) = methods
  
  # All candidate features ranked (except those thrown out for age)
  mutation_dt <- mutation_dt %>%
    arrange(desc(auc))
  
  features_selected <- mutation_dt %>%
    pull(mutation)

  out <- list(features_selected = features_selected,
              select_n = select_n,
              mutation_dt = mutation_dt)
  return(out)
}

# To-do: Write tests