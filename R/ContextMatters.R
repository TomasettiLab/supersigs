# ContextMatters.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 21, 2020
#
# Function for binomial testing using a hierarchical tree structure

# library(dplyr)
# library(rsample)
# library(assertthat)
# library(here)

#' Function for binomial testing using a hierarchical tree structure
#' 
#' Obtain a list of survival mutations by performing a series of nested binomial 
#' tests for selecting and pruning significant mutations
#' 
#' @param muts_df a data frame of mutations
#' @param p_thresh an optional numeric value for the alpha level used in the
#' binomial test (default is `0.05`)
#' @param tot_pseudo an optional numeric value for the pseudo-count value
#' (default is `0`)
#' @param wgs logical value indicating whether sequencing data is 
#' whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE}). 
#' 
#' @import dplyr
#' @import rsample
#' @import assertthat
#' 
#' @export
#' 
#' @return \code{ContextMatters} returns a vector of survival mutations that pass
#' the binomial tree testing
#' 
#' @noRd
#' 
ContextMatters <- function(muts_df, 
                           p_thresh = 0.05,
                           tot_pseudo = 0,
                           wgs = F){
  # Check input
  assert_that(all(setdiff(names(muts_df), "TOTAL_MUTATIONS") == names(muts_formula)),
              msg = "Column names of input data are not correct")
  
  # Use WGS as background probabilities
  if(wgs){
    background_probs = background_probs_wgs
  }
  
  muts_counts <- muts_df %>% colSums() 
  
  # Pseudo-count
  for(feature in names(muts_counts)){
    all_possible_tri["TOTAL_MUTATIONS"] = 3
    muts_counts[feature] = muts_counts[feature] + all_possible_tri[feature]*tot_pseudo/3
  }
  muts_counts <- sapply(muts_counts, round)
  test_all_96 <- T
  bonf_correction <- 150
  
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
# background_probs <- readRDS(here("data", "background_probs.rds"))
# background_probs_wgs <- readRDS(here("data", "background_probs_wgs.rds"))
# muts_children_level3_df <- readRDS(here("data", "muts_children_level3_df.rds"))
# muts_formula <- readRDS(here("data", "muts_formula.rds"))
# h_muts_index <- readRDS(here("data", "h_muts_index.rds"))
# all_possible_tri <- readRDS(here("data", "all_possible_tri.rds"))

# Test function
# signature_caf <- readRDS(here("data", "signature_caf.rds"))
# test_1 = signature_caf[["Data", "ALCOHOL (ESCA)"]]$DataSetFiltered %>%
#   filter(IndVar == 0) %>%
#   transmute_(.dots = muts_formula) %>% # Add counts for every mutation in all_muts
#   mutate(TOTAL_MUTATIONS = select(., 1:6) %>% rowSums())
# 
# ContextMatters(test_1)
