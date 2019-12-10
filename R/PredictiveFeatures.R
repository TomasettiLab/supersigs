# PredictiveFeatures.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Sep 18, 2019
#
# Function for testing and choosing predictive features (formerly part of MyLDAEnvClassifier.R)

library(dplyr)
library(caret)
library(rsample)
library(here)
source(here("code", "MyAuc.R"))

# Helper function to calculate MAD and 25%, 50%, 75% quantiles
get_mad_quantiles <- function(df){
  out <- df %>% apply(MARGIN = 2, function(x)
    c(MAD = mad(x), quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = T))) %>%
    t()
  return(out)
}

# in$train is the subsetted training data frame returned by TransformData
# in$new_partition is the new_partition from GenerateMinSigmaAlgebra
# in$factor is the factor to be tested (e.g. Age, Smoking)
# in$select_n is the number of top features to retain (chosen by cv)
# out$mutation_dt is a data frame of median AUCs and normal p_values (for testing purposes)
# out$features_selected is a vector of predictive features
PredictiveFeatures <- function(train, 
                               new_partition,
                               factor,
                               nmf_out,
                               unsupervised_sig,
                               random_out){
  # Use rates for non-age factors
  if(factor != "AGE"){
    train <- train %>%
      mutate_at(vars(names(new_partition)), funs(./AGE))
  }
  
  # Initialize variables
  train_0 <- train %>% filter(IndVar == 0)
  train_1 <- train %>% filter(IndVar == 1)
  feature_names <- names(new_partition)
  
  # Bootstrap AUC
  boot_dt <- train %>%
    select(c(feature_names, "IndVar")) %>%
    bootstraps(1000, strata = "IndVar")
  
  auc_mat <- sapply(boot_dt$splits, function(x) 
    AllMyAuc(as_tibble(x), DpnVar = "IndVar"))
  
  if(factor == "AGE"){
    auc_mat <- auc_mat  # Only count AUC > 0.5 for age
  } else {
    auc_mat <- abs(auc_mat - 0.5) + 0.5 # Count AUC > 0.5 or < 0.5 for other variables
  }
  
  auc_medians <- apply(auc_mat, 1, median)
  
  # Test which ones are actually predictive (translate AUC to p-values)
  auc_medians <- sapply(auc_medians, function(x) ifelse(is.na(x), 0, x)) # Replace missing AUCs with 0
  n_samples <- nrow(train)
  n1 <- sum(train$IndVar == 1)
  n0 <- sum(train$IndVar == 0)
  
  # Binomial Test
  binom_pvalues <- sapply(auc_medians, function(x)
    qnorm(binom.test(x = floor(n_samples*x), n = n_samples, p = 0.5, alternative = "greater")$p.value,
          lower.tail = F))
  
  # Normal Test
  normal_pvalues <- sapply(auc_medians, function(x)
    qnorm(pnorm(x*n_samples, n_samples/2, sqrt(n_samples)/2, lower.tail = F), lower.tail = F))
  
  # Mann-Whitney-Wilcoxon Test
  # wilcox_pvalues <- qnorm(pwilcox(auc_medians*n0*n1, m = n0, n = n1, lower.tail = F), lower.tail = F)

  # Ad-hoc penalty (this can be improved? e.g. cross-validation)
  feature_lengths <- sapply(new_partition, length)
  penalty <- (log2(96) - log2(feature_lengths))/(2*log2(96)) # penalized z-score, broader feature is better
  
  # Create mutation_dt
  mutation_dt <- as_tibble((get_mad_quantiles(train_0[, feature_names, drop = F]) +
                              get_mad_quantiles(train_1[, feature_names, drop = F]))/2) %>%
    mutate(mutation = feature_names,
           tMad = `50%`/(MAD + 0.1),
           auc = auc_medians,
           binom_pvalues = binom_pvalues, 
           #wilcox_pvalues = wilcox_pvalues,
           normal_pvalues = normal_pvalues,
           normal_penalized = normal_pvalues + (penalty[as.character(mutation)])) %>%
    arrange(desc(normal_pvalues))
  
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
      filter(!(mutation %in% remove_neg_features)) %>%
      arrange(desc(normal_pvalues))
  }
  
  # Inner loop of CV to find the best select_n for each classifier
  max_n = nrow(mutation_dt)
  auc_by_methods_mean = vector("list", max_n)
  methods <- c("LDA", "Logit", "RF")
  for(k in 1:max_n){
    n_iter = 5
    n_fold = 5
    inner_partitions <- sapply(1:n_iter, FUN = function(x) createFolds(y = train$IndVar, k = n_fold))
    inner_aucs <- vector(length = length(inner_partitions), mode = "list")
    
    for(ij in seq_along(inner_partitions)){
      i <- (ij-1) %% n_fold + 1
      j <- floor((ij-1)/n_fold) + 1
      print(paste("Inner CV Iteration:", ij))
      
      features_selected <- mutation_dt %>%
        top_n(n = k, wt = normal_pvalues) %>%
        pull(mutation)
      inner_aucs[[ij]] = MyLDAEnvClassifier(dt = train, 
                         test_ind = inner_partitions[[i, j]],
                         factor,
                         nmf_out,
                         unsupervised_sig,
                         random_out,
                         classifier = methods, 
                         keep_classifier = F,
                         calculate_corr = F,
                         adjusted_formula = F,
                         features_selected = features_selected,
                         select_n = c("LDA" = k, "RF" = k, "Logit" = k))
    }
    
    inner_aucs_mat <- matrix(inner_aucs, nrow = n_fold)
    auc_by_methods <- sapply(methods, function(method) 
      apply(inner_aucs_mat, MARGIN = c(1, 2), function(x){x[[1]]$auc[method]}), simplify = F)
    auc_by_methods_mean[[k]] <- sapply(auc_by_methods, function(x) mean(x, na.rm = T))
  }
  
  select_n = sapply(methods, function(classifier){
    tmp_max = which.max(unlist(lapply(auc_by_methods_mean, function(auc) auc[classifier])))
    if(identical(tmp_max, integer(0))) return(NA) # If classifier fails for all inner folds
    else return(tmp_max)})
  names(select_n) = methods
  
  # All candidate features (except those thrown out for age) ranked
  features_selected <- mutation_dt %>%
    pull(mutation)

  out <- list(mutation_dt = mutation_dt,
              features_selected = features_selected,
              select_n = select_n)
  return(out)
}

# To-do: Write tests