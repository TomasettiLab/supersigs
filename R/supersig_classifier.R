# supersig_classifier.R
# -----------------------------------------------------------------------------
# Author: Bahman Afsari, Albert Kuo
# Date last modified: Feb 25, 2021
#
# Function for classification logistic regression

# library(randomForest)
# library(MASS)
# library(dplyr) 
# library(rsample)
# library(here)
# source(here("code", "MyCor.R"))

# Run logistic regression
run_logit <- function(dt, factor, out, keep_classifier, features_selected,
                      select_n, train_ind, test_ind, test_indvar){
    z <- dt %>%
        select(c(features_selected[seq_len(select_n["Logit"])], "IndVar")) 
    x <- z[train_ind, ]
    newdata <- z[test_ind, ]
    out$auc <- c(out$auc, Logit = NA)
    logit_prediction <- NA
    
    try({
        logit_classifier <- glm(formula = IndVar ~ ., data = x,
                                family = binomial())
        logit_prediction <- predict(logit_classifier, 
                                    newdata = newdata, 
                                    type = "response")
        features_to_summ <- features_selected[seq_len(select_n["Logit"])]
        # Calculate average rates
        if(factor == "AGE"){
            dt_train <- dt %>% slice(train_ind)
            mean_age <- mean(dt$AGE, na.rm = TRUE)
            mean_diffs <- dt %>%
                summarize_at(.vars = features_to_summ,
                             .funs = ~(mean(.)/mean_age))
        } else {
            # Calculate empirical mean differences in rates
            grouped_rates <- dt %>% slice(train_ind) %>% group_by(.data$IndVar)
            grouped_rates <- grouped_rates %>%
                summarize_at(.vars = features_to_summ,
                             .funs = ~(mean(./.data$AGE, na.rm = TRUE)))
            unexposed_rates <- grouped_rates %>% 
                filter(.data$IndVar == FALSE) %>% 
                select(-.data$IndVar)
            exposed_rates <- grouped_rates %>% 
                filter(.data$IndVar == TRUE) %>% 
                select(-.data$IndVar)
            mean_diffs <- exposed_rates - unexposed_rates
        }
        # Save signature for apparent
        if(identical(test_ind, train_ind)){
            names(mean_diffs) <- features_selected[seq_len(select_n["Logit"])]
            
            out$signature <- list(mean_diffs = mean_diffs,
                                  select_n = select_n)
        }
        out$auc["Logit"] <- my_auc(test_indvar, logit_prediction)
        if(keep_classifier) 
            out$classifier$Logit <- logit_classifier
    })
    
    return(out)
}

#' Classification of exposure using signatures
#' 
#' Calculate the AUC using signatures and logistic regression
#' 
#' @param dt a transformed data frame from FeatureSelection
#' @param test_ind an optional vector of indices for the test data
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' @param keep_classifier an optional logical value indicating whether to save
#' the classifier model (default is \code{FALSE})
#' @param features_selected a vector of candidate features ranked by AUC
#' @param select_n the number of top features to retain for each method
#' 
#' @import dplyr
#' @import rsample
#' @importFrom rlang .data
#'
#' @return \code{supersig_classifier} returns a list of several elements:
#' \itemize{
#' \item \code{auc} is a vector of AUCs for each classifier
#' \item \code{signature} is a vector of mean differences or rates
#' (created only when \code{classifier} includes "Logit")
#' \item \code{classifier} is a list of saved models for each input classifier
#' (created only when \code{keep_classifier} is \code{TRUE})
#' }
#' 
#' @noRd
#' 
supersig_classifier <- function(dt, 
                                test_ind = NULL,
                                factor,
                                keep_classifier = FALSE,
                                features_selected,
                                select_n){
        # Split training and test set
        if(is.null(test_ind)){
                train_ind <- seq_len(nrow(dt))
                test_ind <- train_ind
        } else {
                train_ind <- setdiff(seq_len(nrow(dt)), test_ind)
        }
        
        train <- dt[train_ind, ]
        # Initialize output
        test_indvar <- dt[test_ind, "IndVar"]
        out <- list(auc = c(), signature = list(), classifier = list())

        # Transform data
        if(factor != "AGE"){
                dt <- dt %>%
                        mutate_at(.vars = features_selected,
                                  .funs = funs(./.data$AGE))
                
                # Remove observations missing age
                missing_ind <- which(is.na(dt$AGE))
                train_ind <- setdiff(train_ind, missing_ind)
                test_ind <- setdiff(test_ind, missing_ind)
                train <- dt[train_ind, ]
                test_indvar <- dt[test_ind, "IndVar"]
        } 
        
        # Classification
        # If no test data or all training data has the same IndVar
        if(length(test_ind) == 0 || length(unique(train$IndVar)) == 1){
                warning("No valid test data")
                out$auc <- c(NA)
        } else {
            out <- run_logit(dt, factor, out, 
                             keep_classifier, features_selected, 
                             select_n, train_ind, test_ind, test_indvar)
        }
        
        return(out)
}

# No data dependencies

# Test function
# signature_caf <- readRDS(here("data", "signature_caf.rds"))
# factor <- "AGE"
# tissue <- "LUAD"
# ind <- which((signature_caf["Factor",] == factor) & 
# (signature_caf["Tissue",] == tissue))
# dt <- signature_caf[["Data", ind]]$DataSetFiltered %>%
#               filter(TOTAL_MUTATIONS > 0)
# unsupervised_sig = signature_caf[["Unsupervised", ind]]
# age_ind <- which((signature_caf["Factor",] == "AGE") & 
# (signature_caf["Tissue",] == tissue))
# age_sig <- signature_caf[["Unsupervised", age_ind]]
# 
# source(here("code", "FeatureSelection.R"))
# features_out = FeatureSelection(dt = dt, middle_dt = NULL, factor = factor)
# dt = features_out$dt_new
# classifier = c("LDA", "Logit", "RF", "NNLS")
# keep_classifier = F
# features_selected = features_out$features_selected
# select_n = features_out$select_n
# test_ind = NULL
# 
# test_out = supersig_classifier(dt = dt, test_ind = NULL,
#                                factor = factor,
#                                classifier = c("LDA", "Logit", "RF", "NNLS"),
#                                keep_classifier = F,
#                                features_selected = features_selected,
#                                select_n = select_n)

# Create test objects for Shiny app
# test_model = supersig_classifier(dt = dt, test_ind = NULL,
#                                  factor = factor,
#                                  classifier = c("LDA"),
#                                  keep_classifier = T,
#                                  features_selected = features_selected)$classifier
# test_formula = 
# features_out$features_gmsa$new_partition_formula[features_selected]
# saveRDS(list(model = test_model, formula = test_formula), 
# "app/data/test_shiny.rds")
# write.csv(signature_caf[["Data", ind]]$DataSetFiltered[1, ], 
# "app/data/dt.csv", row.names = F)
