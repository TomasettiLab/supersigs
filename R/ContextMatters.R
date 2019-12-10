# ContextMatters.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Aug 6, 2019
#
# Function for binomial testing using a hierarchical tree structure

library(dplyr)
library(rsample)
library(assertthat)
library(here)

# in$muts_df is the dataset of mutations
# in$p_thresh is the alpha level for testing
# in$min_median and in$min_samples is the minimum threshold to test 96 instead of 6 mutations (deprecated)
# in$pseudo_counts is a toggle for the bootstrap + pseudocount methodology
# output is a vector of survival mutations that passed the tree binomial testing
ContextMatters <- function(muts_df, p_thresh = 0.05, 
                           min_median = 20, min_samples = 60,
                           use_pseudo = T){
  # Check input
  assert_that(all(setdiff(names(muts_df), "TOTAL_MUTATIONS") == names(muts_formula)),
              msg = "Column names of input data are not correct")
  
  if(use_pseudo){
    bootstrapped_df <- muts_df %>% select(c("TOTAL_MUTATIONS", names(muts_formula))) %>%
      bootstraps(times = 100)
    
    # Bootstrap
    muts_counts <- bootstrapped_df$splits %>%
      sapply(FUN = function(x){
        z <- as_tibble(x) %>% colSums()}) %>%
      apply(MARGIN = 1, FUN = median, na.rm = T)
   
    # Pseudo-count
    tot_pseudo <- 1000
    for(feature in names(muts_counts)){
      all_possible_tri["TOTAL_MUTATIONS"] = 3
      muts_counts[feature] = muts_counts[feature] + all_possible_tri[feature]*tot_pseudo/3
    }
    muts_counts <- sapply(muts_counts, round)
    test_all_96 <- T
    bonf_correction <- 150
  } else {
    muts_counts <- muts_df %>% colSums() 
    test_all_96 <- (median(muts_df$TOTAL_MUTATIONS, na.rm=T) > min_median) && 
      (nrow(muts_df) >= min_samples)
    bonf_correction <- 150
    
    # Limit feature space to the 6 central mutations if minimums not met
    if(!test_all_96){
      background_probs <- background_probs %>%
        filter(feature %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
      bonf_correction <- 6
    }
  }
  
  # Create "tree"
  tree <- tibble(feature = names(h_muts_index),
                 n_children = h_muts_index) %>%
    right_join(. , background_probs, by = "feature") %>%
    mutate(q = muts_counts[feature],
           size = muts_counts[parent_name],
           prop = q/size)
  
  # Binomial test for all features
  tree <- tree %>%
    mutate(p_value = pbinom(q, size, prob, lower.tail = F), 
           #p_value_bonf = p.adjust(p_value, method = "bonferroni", n = nrow(tree)),  # Supplement says 151 features, but there are 486 = nrow(tree) tests
           p_value_bonf = p_value*bonf_correction,
           sig = p_value_bonf < p_thresh) 
  
  # First level of significant features
  survival_mutations <- tree %>%
    filter(n_children == 16) %>%
    filter(sig) %>% 
    pull(feature)
  
  # Return here if only testing for 6 central mutations
  if(!test_all_96){
    return(survival_mutations)
  }
  
  # Keep tests only for significant parents
  for(i in c(4, 1)){
    tree_tmp <- tree %>%
      filter(n_children == i) %>%
      mutate(sig_2 = ifelse(parent_name %in% c(survival_mutations, "TOTAL_MUTATIONS"), sig, T))
    
    survival_mutations_tmp <- tree_tmp %>% 
      filter(n_children == i) %>%
      group_by(feature) %>%
      summarize(sig_2 = all(sig_2)) %>%
      filter(sig_2) %>%
      pull(feature)
    
    survival_mutations <- c(survival_mutations, unique(survival_mutations_tmp))
  }
  
  survival_mutations_cache <- tree %>%
    filter(feature %in% survival_mutations)
  
  # Keep significance only if still significant
  # when all significant children contributions are removed
  for(level in 2:1){
    # Take union of children and grandchildren
    union_children_vec = survival_mutations_cache %>%
      filter(n_children %in% if(level == 1) c(4, 1) else c(1)) %>%
      left_join(. , muts_children_level3_df, by = "feature") %>%
      pull(child_name) %>%
      unique()
    
    # All survival children and grandchildren as level 3 mutations
    survival_children_level3 <- tree %>%
      filter(feature %in% union_children_vec) %>%
      select(feature, parent_name, prob, q) %>%
      rename(child_name = feature,
             feature = parent_name,
             prob_child = prob,
             q_child = q)
    
    # Compute children of feature
    tree_pruned <- tree %>%
      filter(feature %in% survival_mutations &
               n_children == ifelse(level == 1, 16, 4)) %>%
      select(feature, n_children, parent_name, prob, q, size) %>%
      left_join(. , muts_children_level3_df, by = "feature") %>%
      inner_join(. , survival_children_level3, by = c("child_name", "parent_name" = "feature")) %>%
      filter(parent_name %in% c(survival_mutations, "TOTAL_MUTATIONS")) %>%
      group_by(feature, parent_name, prob, q, size) %>%
      summarize(prob_child = sum(prob_child),
                q_child = sum(q_child)) %>%
      ungroup()
    
    # Compute children of parent
    tree_pruned_parent <- tree_pruned %>%
      distinct(parent_name) %>%
      rename(feature = parent_name) %>%
      left_join(. , muts_children_level3_df, by = "feature") %>%
      inner_join(. , survival_children_level3, by = c("child_name", "feature")) %>%
      group_by(feature) %>%
      summarize(prob_parent_child = sum(prob_child),
                q_parent_child = sum(q_child)) %>%
      ungroup()
    
    # Join and test
    tree_pruned <- tree_pruned %>%
      full_join(. , tree_pruned_parent, by = c("parent_name" = "feature")) %>%
      mutate(prob_pruned = (prob - prob_child)/(1 - prob_parent_child),
             q_pruned = q - q_child,
             size_pruned = size - q_parent_child) %>%
      mutate(p_value = pbinom(q_pruned, size_pruned, prob_pruned, lower.tail = F), 
             p_value_bonf = p_value*bonf_correction,
             sig = p_value_bonf < p_thresh)
    
    pruned_mutations <- tree_pruned %>%
      group_by(feature) %>%
      summarize(sig = all(sig)) %>%
      filter(!sig) %>%
      pull(feature)
    
    survival_mutations_cache <- survival_mutations_cache %>%
      filter(!(feature %in% pruned_mutations))
    
    survival_mutations <- unique(survival_mutations_cache$feature)
  }
  
  if(length(survival_mutations) == 0){
    return("TOTAL_MUTATIONS")
  } else {
    return(survival_mutations)
  }
}

# Load data dependencies
background_probs <- readRDS(here("data", "background_probs.rds"))
muts_children_level3_df <- readRDS(here("data", "muts_children_level3_df.rds"))
muts_formula <- readRDS(here("data", "muts_formula.rds"))
h_muts_index <- readRDS(here("data", "h_muts_index.rds"))
all_possible_tri <- readRDS(here("data", "all_possible_tri.rds"))

# Test function
# signature_caf <- readRDS(here("data", "signature_caf.rds"))
# test_1 = signature_caf[["Data", "ALCOHOL (ESCA)"]]$DataSetFiltered %>%
#   filter(IndVar == 0) %>%
#   transmute_(.dots = muts_formula) %>% # Add counts for every mutation in all_muts
#   mutate(TOTAL_MUTATIONS = select(., 1:6) %>% rowSums())
# 
# ContextMatters(test_1)
