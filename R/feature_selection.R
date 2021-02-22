# feature_selection.R
# -----------------------------------------------------------------------------
# Author: Bahman Afsari, Albert Kuo
# Date last modified: Feb 22, 2021
#
# Function for selecting features

# library(dplyr)
# library(assertthat)
# library(caret)
# library(here)
# source(here("R", "context_matters.R"))
# source(here("R", "generate_min_sigma_algebra.R"))
# source(here("R", "transform_data.R"))

# Add up counts for every mutation (helper function)
add_counts <- function(dt){
    dt <- dt %>%
        transmute_(.dots = muts_formula) %>%
        mutate(TOTAL_MUTATIONS = 
                         select(., c("C>A", "C>G", "C>T", 
                                     "T>A", "T>C", "T>G")) %>%
                         rowSums())
    
    return(dt)
}

# Apply context matters to every inner partition (helper function)
feature_engineering <- function(train, factor, inner_partitions, n_fold, wgs){
    message(paste("Begin feature engineering..."))
    new_partition_ls <- vector(length = length(inner_partitions), 
                               mode = "list")
    for(ij in seq_along(inner_partitions)){
        i <- (ij-1) %% n_fold + 1
        j <- floor((ij-1)/n_fold) + 1
        
        test_inner_ind <- inner_partitions[[i, j]]
        train_inner_ind <- setdiff(seq_len(nrow(train)), test_inner_ind)
        train_inner <- train[train_inner_ind, ]

        # Separate into exposed (train_1) and unexposed (train_0)
        if(factor != "AGE"){ 
            train_0 <- train_inner %>% filter(.data$IndVar == 0)
            train_1 <- train_inner %>% filter(.data$IndVar == 1)
        } else {
            train_0 <- train_inner 
            train_1 <- train_0
        }
        
        # Add up counts for every mutation
        train_0 <- add_counts(train_0)
        train_1 <- add_counts(train_1)

        # Test for significant features
        if(nrow(train_0) == 0 || nrow(train_1) == 0) next
        features_context_0 <- context_matters(train_0, wgs = wgs)
        features_context_1 <- context_matters(train_1, wgs = wgs)
        
        input_ls <- list(var0 = features_context_0, var1 = features_context_1)
        features_gmsa <- generate_min_sigma_algebra(input_ls)
        new_partition_ls[[ij]] <- features_gmsa
    }
    
    # Repartition features
    new_partition_combined <- Reduce(function(a, b){
        generate_min_sigma_algebra(list(var0 = a, var1 = b),
                                   partitioned_features = TRUE)},
        new_partition_ls)
    
    return(new_partition_combined)
}

# Rank features by AUC (helper function)
rank_tri_features <- function(dt_new, inner_partitions,
                              new_partition, factor){
    rank_tri_ls <- vector(length = length(inner_partitions), mode = "list")
    for(ij in seq_along(inner_partitions)){
        feature_names <- names(new_partition)
        auc_mat <- all_my_auc(dt_new %>% select(c(feature_names, "IndVar")), 
                                                    IndVar = "IndVar")
        if(factor == "AGE"){
            auc_medians <- auc_mat
        } else {
            auc_medians <- ifelse(auc_mat > 0.5, auc_mat, 1 - auc_mat)
        }
        rank_mutations <- tibble(mutation = feature_names,
                                 auc = auc_medians) %>%
            arrange(desc(.data$auc))
        
        # Calculate rank of trinucleotide equivalents
        feature_to_trinuc <- tibble(feature = rep(feature_names,
                                                  vapply(new_partition, 
                                                                 length,
                                                                 numeric(1))),
                                    trinucleotide = unlist(new_partition))
        rank_trinucleotides <- rank_mutations %>%
            left_join(., feature_to_trinuc,
                                by = c("mutation" = "feature")) %>%
            mutate(rank = seq_len(n())) %>%
            group_by(.data$mutation) %>%
            mutate(mean_rank = mean(.data$rank)) %>%
            ungroup()
        
        # Save trinucleotide rankings
        rank_tri_ls[[ij]] <- rank_trinucleotides %>%
            select(.data$trinucleotide, .data$mean_rank) %>%
            arrange(.data$trinucleotide) %>%
            dplyr::rename(!!paste0("rank_", ij) := .data$mean_rank,
                          !!paste0("trinucleotide_", ij) := .data$trinucleotide)
    }
    
    return(rank_tri_ls)
}

# Get repartition of features (helper function)
repartition_with_rank <- function(rank_tri_ls){
    # Take median of rankings for trinucleotides over inner folds
    rank_tri <- bind_cols(rank_tri_ls)
    rank_tri <- rank_tri %>%
        mutate(rank = apply(select(., starts_with("rank")), 1, median)) %>%
        select(.data$trinucleotide_1, .data$rank) %>%
        dplyr::rename(trinucleotide = .data$trinucleotide_1) %>%
        arrange(.data$rank)
    
    # Create new partition of features, grouping trinucleotides of same rank
    # into one feature
    unique_ranks <- unique(rank_tri$rank)
    features_new <- lapply(unique_ranks, function(r) rank_tri %>%
                             filter(.data$rank == r) %>%
                             pull(.data$trinucleotide)) 
    names(features_new) <- paste0("X", seq_along(features_new))
    new_partition <- features_new
    features_gmsa <- list(new_partition = features_new)
    features_selected <- names(new_partition)
    
    return(list(features_gmsa, new_partition, features_selected))
}

# Get n star (helper function)
n_star_features <- function(train, factor, inner_partitions,
                            n_fold, features_selected){
    n_star_ls <- vector(length = length(inner_partitions), mode = "list")
    
    n <- length(features_selected)
    message(paste("Begin cross-validated selection over", n, "features and", 
                  length(inner_partitions), "inner folds..."))
    
    for(ij in seq_along(inner_partitions)){
        message(paste("...testing inner fold", ij))
        i <- (ij-1) %% n_fold + 1
        j <- floor((ij-1)/n_fold) + 1
        
        test_inner_ind <- inner_partitions[[i, j]]
        train_inner_ind <- setdiff(seq_len(nrow(train)), test_inner_ind)
        
        # Find n*
        auc_ls <- list()
        methods <- c("Logit")
        
        for(k in seq_len(n)){
            auc_ls[[k]] <- supersig_classifier(dt = train, 
                                               test_ind = test_inner_ind,
                                               factor,
                                               keep_classifier = FALSE,
                                               features_selected,
                                               select_n = c("Logit" = k))$auc
        }
        n_star <- which.max(vapply(auc_ls, function(x) x["Logit"],
                                   FUN.VALUE = numeric(1)))
        n_star <- ifelse(identical(n_star, integer(0)), NA, n_star)
        
        # Save n_star for each method
        n_star_ls[[ij]] <- tibble(!!paste0("methods_", ij) := methods,
                                  !!paste0("n_star_", ij) := n_star)
    }
    
    return(n_star_ls)
}

# Get median n star (helper function)
median_n_star <- function(n_star_ls){
    n_star <- bind_cols(n_star_ls)
    n_star <- n_star %>%
        mutate(n_star = apply(select(., starts_with("n_star")), 1, 
                              median, na.rm = TRUE),
               n_star = round(.data$n_star)) %>%
        select(.data$methods_1, .data$n_star) %>%
        dplyr::rename(methods = .data$methods_1)
    
    select_n <- n_star$n_star
    names(select_n) <- n_star$methods
    
    return(select_n)
}

# Adjust select n by removing features that are not > 0.6 AUC
select_n_adjust <- function(train, features_selected, select_n){
    auc_mat <- all_my_auc(train %>% select(c(features_selected, "IndVar")), 
                          IndVar = "IndVar")
    max_n <- sum(auc_mat > 0.6)
    select_n <- vapply(select_n, function(x) min(max_n, x),
                       FUN.VALUE = numeric(1))
}

#' Function for feature selection
#' 
#' Perform feature selection given a dataset of mutations by calling 
#' a series of other functions to find significant and predictive features
#' 
#' @param dt a data frame of mutations
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' #' @param wgs logical value indicating whether sequencing data is 
#' whole-genome (wgs = \code{TRUE}) or whole-exome (wgs = \code{FALSE}). 
#' 
#' @import dplyr
#' @import assertthat
#' @import caret
#' @importFrom rlang .data
#' 
#' @return \code{feature_selection} returns a list of several elements:
#' \itemize{
#' \item \code{features_context_0} is a vector of survival mutations
#' for the unexposed group
#' \item \code{features_context_1} is a vector of survival mutations
#' for the exposed group
#' \item \code{features_selected} is a vector of candidate 
#' features ranked by AUC
#' \item \code{new_partition} is a list of the partitioned candidate features,
#' where each feature is represented by a vector of fundamental mutations
#' \item \code{select_n} is the number of top features to retain for each method
#' \item \code{mutation_dt} is a sorted data frame of candidate features 
#' and their AUC
#' \item \code{dt_new} is the transformed data from \code{transform_data}
#' }
#' 
#' @noRd
#' 
feature_selection <- function(dt,
                              factor,
                              wgs = FALSE){
        train <- dt
        n_iter <- 5
        n_fold <- 3
        inner_partitions <- vapply(seq_len(n_iter), 
                                   FUN = function(x){
                                     createFolds(y = train$IndVar,
                                                 k = n_fold)},
                                   FUN.VALUE = vector("list", length = n_fold))
        
        # Feature engineering
        new_partition <- feature_engineering(train, factor, inner_partitions, 
                                             n_fold, wgs)
        
        # Transform data
        dt_new <- transform_data(train, new_partition)
        
        # Trinucleotide ranking in inner folds
        rank_tri_ls <- rank_tri_features(dt_new, inner_partitions,
                                         new_partition, factor)
        
        # Get new partition after grouping trinucleotides of same rank
        repartition <- repartition_with_rank(rank_tri_ls)
        features_gmsa <- repartition[[1]]
        new_partition <- repartition[[2]]
        features_selected <- repartition[[3]]
        
        # Transform data 
        dt_new <- transform_data(dt, new_partition)
        train <- dt_new
        
        # Find n*
        n_star_ls <- n_star_features(train, factor, inner_partitions, 
                                     n_fold, features_selected)
        select_n <- median_n_star(n_star_ls) # take medians over inner folds
        select_n <- select_n_adjust(train, features_selected, select_n)
        
        # Return output
        assert_that(length(features_selected) >= 1, 
                    msg = "No predictive features found")
        
        out <- list(features_gmsa = features_gmsa, 
                    features_selected = features_selected,
                    new_partition = new_partition,
                    dt_new = dt_new,
                    select_n = select_n)
        return(out)
}

# Load data dependencies
# muts_formula <- readRDS(here("data", "muts_formula.rds"))

# Test function
# signature_caf <- readRDS(here("super_sigs", "data", "signature_caf.rds"))
# factor <- "AGE"
# tissue <- "UCEC"
# ind <- which((signature_caf["Factor",] == factor) &
# (signature_caf["Tissue",] == tissue))
# dt <- signature_caf[["Data", ind]]$DataSetFiltered %>%
#               filter(TOTAL_MUTATIONS > 0)
# 
# test_ind = NULL
# middle_dt = NULL
# 
# if(factor == "AGE"){
#               middle_dt <- signature_caf[["Data", ind]]$DataSetFilteredKeepMiddle %>%
#                       filter(!(SAMPLE %in% dt$SAMPLE)) %>%
#                       mutate(IndVar = NA)
# }
# 
# test_out = feature_selection(dt = dt, middle_dt = middle_dt, factor = factor)