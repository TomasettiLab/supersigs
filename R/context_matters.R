# context_matters.R
# -----------------------------------------------------------------------------
# Author: Bahman Afsari, Albert Kuo
# Date last modified: Feb 18, 2021
#
# Function for binomial testing using a hierarchical tree structure

# library(dplyr)
# library(rsample)
# library(assertthat)
# library(here)


# Binomial test for all features in tree (helper function)
binom_test_tree <- function(tree){
    p_thresh <- 0.05
    bonf_correction <- 150
    
    tree <- tree %>%
        mutate(p_value = pbinom(.data$q, .data$size, .data$prob,
                                lower.tail = FALSE), 
               p_value = ifelse(.data$size == 0, 1, .data$p_value),
               p_value_bonf = .data$p_value*bonf_correction,
               sig = .data$p_value_bonf < p_thresh) 
    
    return(tree)
}

# Create tree (helper function)
create_tree <- function(muts_counts, background_probs){
    # Create "tree"
    tree <- tibble(feature = names(h_muts_index),
                   n_children = h_muts_index) %>%
        right_join(. , background_probs, by = "feature") %>%
        mutate(q = muts_counts[feature],
               size = muts_counts[.data$parent_name],
               prop = .data$q/.data$size)
    
    # Binomial test for all features
    tree <- binom_test_tree(tree)
    
    # First level of significant features
    survival_mutations <- tree %>%
        filter(.data$n_children == 16) %>%
        filter(.data$sig) %>% 
        pull(feature)
    
    # Keep tests only for significant parents
    for(i in c(4, 1)){
        tree_tmp <- tree %>%
            filter(.data$n_children == i) %>%
            mutate(sig_2 = ifelse(.data$parent_name %in% 
                                    c(survival_mutations, "TOTAL_MUTATIONS"), 
                                  .data$sig, TRUE))
        
        survival_mutations_tmp <- tree_tmp %>% 
            filter(.data$n_children == i) %>%
            group_by(.data$feature) %>%
            summarize(sig_2 = all(.data$sig_2)) %>%
            filter(.data$sig_2) %>%
            pull(.data$feature)
        
        survival_mutations <- c(survival_mutations, 
                                unique(survival_mutations_tmp))
    }
    
    survival_mutations_cache <- tree %>%
        filter(feature %in% survival_mutations)
    
    return(list(tree, survival_mutations_cache, survival_mutations))
}

# Compute children of feature (helper function)
children_of_feature <- function(tree, survival_mutations, 
                                survival_children_level3, level){
    tree_pruned <- tree %>%
        filter(.data$feature %in% survival_mutations &
                 .data$n_children == ifelse(level == 1, 16, 4)) %>%
        select(.data$feature, .data$n_children, .data$parent_name, 
               .data$prob, .data$q, .data$size) %>%
        left_join(. , muts_children_level3_df, by = "feature") %>%
        inner_join(. , survival_children_level3, 
                   by = c("child_name", "parent_name" = "feature")) %>%
        filter(.data$parent_name %in% 
                 c(survival_mutations, "TOTAL_MUTATIONS")) %>%
        group_by(.data$feature, .data$parent_name, .data$prob, 
                 .data$q, .data$size) %>%
        summarize(prob_child = sum(.data$prob_child),
                  q_child = sum(.data$q_child)) %>%
        ungroup()
    
    return(tree_pruned)
}

# Compute children of parent (helper function)
children_of_parent <- function(tree_pruned, survival_children_level3){
    tree_pruned_parent <- tree_pruned %>%
        distinct(.data$parent_name) %>%
        rename(feature = .data$parent_name) %>%
        left_join(. , muts_children_level3_df, by = "feature") %>%
        inner_join(. , survival_children_level3,
                   by = c("child_name", "feature")) %>%
        group_by(.data$feature) %>%
        summarize(prob_parent_child = sum(.data$prob_child),
                  q_parent_child = sum(.data$q_child)) %>%
        ungroup()
    
    return(tree_pruned_parent)
}

# Join and test pruned tree (helper function)
test_pruned_tree <- function(tree_pruned, tree_pruned_parent){
    p_thresh <- 0.05
    bonf_correction <- 150
    
    tree_pruned <- tree_pruned %>%
        full_join(. , tree_pruned_parent,
                  by = c("parent_name" = "feature")) %>%
        mutate(prob_pruned = (.data$prob - .data$prob_child)/
                         (1 - .data$prob_parent_child),
               q_pruned = .data$q - .data$q_child,
               size_pruned = .data$size - .data$q_parent_child) %>%
        mutate(p_value = pbinom(.data$q_pruned, .data$size_pruned, 
                                .data$prob_pruned, lower.tail = FALSE), 
               p_value = ifelse(.data$size_pruned == 0, 1, .data$p_value),
               p_value_bonf = .data$p_value*bonf_correction,
               sig = .data$p_value_bonf < p_thresh)
    return(tree_pruned)
}


# Retest tree with significant children contributions removed (helper function)
retest_tree <- function(tree, survival_mutations_cache, survival_mutations){
    p_thresh <- 0.05
    bonf_correction <- 150
    
    for(level in 2:1){
        # Take union of children and grandchildren
        union_children_vec <- survival_mutations_cache %>%
            filter(.data$n_children %in% if(level == 1) c(4, 1) else c(1)) %>%
            left_join(. , muts_children_level3_df, by = "feature") %>%
            pull(.data$child_name) %>%
            unique()
        
        # All survival children and grandchildren as level 3 mutations
        survival_children_level3 <- tree %>%
            filter(.data$feature %in% union_children_vec) %>%
            select(.data$feature, .data$parent_name, .data$prob, .data$q) %>%
            rename(child_name = .data$feature,
                   feature = .data$parent_name,
                   prob_child = .data$prob,
                   q_child = .data$q)
        
        # Compute children of feature
        tree_pruned <- children_of_feature(tree, survival_mutations, 
                                           survival_children_level3, level)
        
        # Compute children of parent
        tree_pruned_parent <- children_of_parent(tree_pruned, 
                                                 survival_children_level3)

        # Join and test
        tree_pruned <- test_pruned_tree(tree_pruned, tree_pruned_parent)
        
        pruned_mutations <- tree_pruned %>%
            group_by(.data$feature) %>%
            summarize(sig = all(.data$sig)) %>%
            filter(!.data$sig) %>%
            pull(.data$feature)
        
        survival_mutations_cache <- survival_mutations_cache %>%
            filter(!(.data$feature %in% pruned_mutations))
        
        survival_mutations <- unique(survival_mutations_cache$feature)
    }
    
    return(survival_mutations)
}

#' Function for binomial testing using a hierarchical tree structure
#' 
#' Obtain a list of survival mutations by performing a series of nested 
#' binomial tests for selecting and pruning significant mutations
#' 
#' @param muts_df a data frame of mutations
#' @param tot_pseudo an optional numeric value for the pseudo-count value
#' (default is `0`)
#' @param wgs logical value indicating whether sequencing data is 
#' whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE}). 
#' 
#' @import dplyr
#' @import assertthat
#' @importFrom rlang .data
#' 
#' @return \code{context_matters} returns a vector of survival mutations that
#' pass the binomial tree testing
#' 
#' @noRd
#' 
context_matters <- function(muts_df, 
                            tot_pseudo = 0, 
                            wgs = FALSE){
    # Check input
    assert_that(all(setdiff(names(muts_df), "TOTAL_MUTATIONS") == 
                      names(muts_formula)),
                      msg = "Column names of input data are not correct 
                            in context_matters")
    
    # Use WGS as background probabilities
    if(wgs) background_probs <- background_probs_wgs
    
    muts_counts <- muts_df %>% colSums()
    # Pseudo-count
    for(feature in names(muts_counts)){
        all_possible_tri["TOTAL_MUTATIONS"] <- 3
        muts_counts[feature] <- muts_counts[feature] + 
            all_possible_tri[feature]*tot_pseudo/3
    }
    muts_counts <- vapply(muts_counts, round, FUN.VALUE = numeric(1))

    # Create tree
    tree_ls <- create_tree(muts_counts, background_probs)
    tree <- tree_ls[[1]]
    survival_mutations_cache <- tree_ls[[2]]
    survival_mutations <- tree_ls[[3]]
    
    # Retest tree after removing significant children
    survival_mutations <- retest_tree(tree, survival_mutations_cache, 
                                      survival_mutations)

    if(length(survival_mutations) == 0){
        return("TOTAL_MUTATIONS")
    } else {
        return(survival_mutations)
    }
}

# Load data dependencies
# background_probs <- readRDS(here("data", "background_probs.rds"))
# background_probs_wgs <- readRDS(here("data", "background_probs_wgs.rds"))
# muts_children_level3_df <- 
# readRDS(here("data", "muts_children_level3_df.rds"))
# muts_formula <- readRDS(here("data", "muts_formula.rds"))
# h_muts_index <- readRDS(here("data", "h_muts_index.rds"))
# all_possible_tri <- readRDS(here("data", "all_possible_tri.rds"))

# Test function
# signature_caf <- readRDS(here("data", "signature_caf.rds"))
# test_1 = signature_caf[["Data", "ALCOHOL (ESCA)"]]$DataSetFiltered %>%
#       filter(IndVar == 0) %>%
#       transmute_(.dots = muts_formula) %>% 
#       mutate(TOTAL_MUTATIONS = select(., 1:6) %>% rowSums())
# 
# context_matters(test_1)
