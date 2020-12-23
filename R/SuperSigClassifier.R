# SuperSigClassifier.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019
#
# Function for classification and correlation using LDA, logistic, or random forest

# library(randomForest)
# library(MASS)
# library(dplyr) # Note that you must load dplyr after MASS to preserve select function
# library(rsample)
# library(here)
# source(here("code", "MyCor.R"))

#' Classification of exposure using signatures
#' 
#' Calculate the AUC using signatures with the option of three 
#' different classifiers
#' 
#' @param dt a transformed data frame from FeatureSelection
#' @param test_ind an optional vector of indices for the test data
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' @param classifier a vector of the classifier method(s) to use for prediction
#' (options are "LDA", "Logit", and "RF")
#' @param keep_classifier an optional logical value indicating whether to save
#' the classifier model (default is \code{FALSE})
#' @param adjusted_formula an optional logical value indicating whether to
#' use the adjusted formula for non-age factors (default is \code{FALSE})
#' @param features_selected a vector of candidate features ranked by AUC
#' @param select_n the number of top features to retain for each method
#' 
#' @import dplyr
#' @import rsample
#' @importFrom randomForest randomForest
#' @importFrom MASS lda
#'
#' @return \code{SuperSigClassifier} returns a list of several elements:
#' \itemize{
#' \item \code{auc} is a vector of AUCs for each classifier
#' \item \code{signature} is a vector of mean differences or rates
#' (created only when \code{classifier} includes "Logit")
#' \item \code{classifier} is a list of saved models for each input classifier
#' (created only when \code{keep_classifier} is \code{TRUE})
#' }
#' 
SuperSigClassifier <- function(dt, test_ind = NULL,
                               factor,
                               classifier, # Options are "LDA", "Logit", "RF"
                               keep_classifier = F,
                               adjusted_formula = F,
                               features_selected,
                               select_n){
  # Split training and test set
  if(is.null(test_ind)){
    train_ind <- 1:nrow(dt)
    test_ind <- train_ind
  } else {
    train_ind <- setdiff(1:nrow(dt), test_ind)
  }
  
  train <- dt[train_ind, ]
  
  # Do not perform LDA if select_n["LDA"] is missing
  if(is.na(select_n["LDA"])){
    classifier = setdiff(classifier, "LDA")
  }
  
  # Initialize output
  test_indvar <- dt[test_ind, "IndVar"]
  out <- list(auc = c(),
              signature = list(),
              classifier = list())

  # Transform data
  if(factor != "AGE"){
    if(adjusted_formula){
      # Calculate grouped median rates
      grouped_rates <- train %>%
        group_by(IndVar) %>%
        summarize_at(.vars = features_selected,
                     .funs = funs(median(., na.rm=T)/median(AGE, na.rm=T)))
      
      unexposed_rates <- grouped_rates %>% filter(IndVar == F) %>% select(-IndVar)
      exposed_rates <- grouped_rates %>% filter(IndVar == T) %>% select(-IndVar)
      
      # Remove unexposed median rate (i.e. aging rate) from training and test data
      remove_age_formula <- colnames(unexposed_rates) %>%
        sapply(FUN = function(x) paste0("`", x, "`", "-AGE*", unexposed_rates[x]))
      
      dt <- dt %>%
        mutate_(.dots = remove_age_formula) %>%
        mutate_at(.vars = features_selected,
                  .funs = funs(./TOTAL_MUTATIONS))
    } else {
      dt <- dt %>%
        mutate_at(.vars = features_selected,
                  .funs = funs(./AGE))
    }
    
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
    warning("No test data")
    
    # Return NA to AUC
    auc <- sapply(classifier, function(x) NA)
    out$auc <- c(out$auc, auc)
    
    # Return NA to corr
    if(factor %in% c("AGE", "SMOKING")){
      corr <- sapply(classifier, function(x) NA)
      out$corr <- c(out$corr, corr)
    }
  } else {
    # Linear Discriminant Analysis
    if("LDA" %in% classifier){
      if(factor == "AGE"){
        # Take log2 for counts for age for better normal distribution
        # dt <- dt %>%
        #   mutate_at(.vars = features_selected[1:select_n["LDA"]],
        #             .funs = funs(log2(. + 1)))
        
        train <- dt[train_ind, ]
      }
      
      x <- train[, features_selected[1:select_n["LDA"]], drop = F] %>% as.matrix()
      grouping <- train$IndVar
      out$auc <- c(out$auc, LDA = NA)
      lda_prediction <- NA
      lda_rates <- NA
      try({
        lda_classifier <- lda(x, grouping, prior = c(0.5, 0.5))
        newdata <- dt[test_ind, features_selected[1:select_n["LDA"]], drop = F] %>% as.matrix()
        lda_prediction <- predict(lda_classifier, newdata = newdata)$x
        
        out$auc["LDA"] = MyAuc(test_indvar, lda_prediction)
        if(keep_classifier) 
          out$classifier$LDA <- lda_classifier
        
        # Calculate mean difference in rates or counts
        lda_rates <- lda_classifier$means %>% as.data.frame()
        lda_rates <- lda_rates %>% mutate(IndVar = as.logical(rownames(lda_rates)))
        unexposed_lda_rates <- lda_rates %>% filter(IndVar == F) %>% select(-IndVar)
        exposed_lda_rates <- lda_rates %>% filter(IndVar == T) %>% select(-IndVar)
        lda_rates <- exposed_lda_rates - unexposed_lda_rates
      })
    }
    
    # Random Forest
    if("RF" %in% classifier){
      x <- train[, features_selected[1:select_n["RF"]], drop = F] %>% as.matrix()
      y_factor <- as.factor(train$IndVar)
      out$auc <- c(out$auc, RF = NA)
      rf_prediction <- NA
      try({
        rf_classifier <- randomForest(x, y = y_factor,
                                      mtry = min(ncol(x), round(log2(ncol(dt)))),
                                      maxnodes = 2, classwt = c(0.5, 0.5))
        
        newdata <- dt[test_ind, features_selected[1:select_n["RF"]], drop = F] %>% as.matrix()
        rf_prediction <- predict(rf_classifier, newdata = newdata,
                                 type = "prob")[, levels(y_factor)[2]]
        
        out$auc["RF"] <- MyAuc(test_indvar, rf_prediction)
        if(keep_classifier) 
          out$classifier$RF <- rf_classifier
      })
    }
    
    # Logistic Regression
    if("Logit" %in% classifier){
      z <- dt %>%
        select(c(features_selected[1:select_n["Logit"]], "IndVar")) # %>%
        # rename_at(vars(features_selected[1:select_n["Logit"]]), function(x) paste0("X", 1:length(features_selected[1:select_n["Logit"]])))
      
      x <- z[train_ind, ]
      newdata <- z[test_ind, ]
      out$auc <- c(out$auc, Logit = NA)
      logit_prediction <- NA
      try({
        logit_classifier <- glm(formula = IndVar ~ ., data = x, family = binomial())
        logit_prediction <- predict(logit_classifier, newdata = newdata, type = "response")
        logit_betas <- coef(logit_classifier) # beta vector
        logit_betas <- coef(logit_classifier)[-1] # beta vector except beta_0
        
        # Calculate empirical mean differences
        if(factor == "AGE"){
          grouped_rates <- dt %>%
            slice(train_ind) %>%
            group_by(IndVar) %>%
            summarize_at(.vars = features_selected[1:select_n["Logit"]],
                         .funs = funs(mean(., na.rm = T)))
        } else {
          grouped_rates <- dt %>%
            slice(train_ind) %>%
            group_by(IndVar) %>%
            summarize_at(.vars = features_selected[1:select_n["Logit"]],
                         .funs = funs(mean(./AGE, na.rm = T)))
        }
        
        unexposed_rates <- grouped_rates %>% filter(IndVar == F) %>% select(-IndVar)
        exposed_rates <- grouped_rates %>% filter(IndVar == T) %>% select(-IndVar)
        
        # Calculate mean_diff = difference in counts or rates between exposed and unexposed -> signature representation
        mean_diffs <- exposed_rates - unexposed_rates
        
        # Save signature (empirical mean rates and beta coefficients) and select_n for apparent
        if(identical(test_ind, train_ind)){
          # Logit beta coefficients
          names(logit_betas) <- features_selected[1:select_n["Logit"]]
          names(mean_diffs) <- features_selected[1:select_n["Logit"]]
          
          out$signature <- list(mean_diffs = mean_diffs, logit_betas = logit_betas, select_n = select_n)
        }
        
        out$auc["Logit"] = MyAuc(test_indvar, logit_prediction)
        if(keep_classifier) 
          out$classifier$Logit <- logit_classifier
      })
    }
  }
  
  return(out)
}


# No data dependencies

# Test function
# signature_caf <- readRDS(here("data", "signature_caf.rds"))
# factor <- "AGE"
# tissue <- "LUAD"
# ind <- which((signature_caf["Factor",] == factor) & (signature_caf["Tissue",] == tissue))
# dt <- signature_caf[["Data", ind]]$DataSetFiltered %>%
#   filter(TOTAL_MUTATIONS > 0)
# unsupervised_sig = signature_caf[["Unsupervised", ind]]
# age_ind <- which((signature_caf["Factor",] == "AGE") & (signature_caf["Tissue",] == tissue))
# age_sig <- signature_caf[["Unsupervised", age_ind]]
# 
# source(here("code", "FeatureSelection.R"))
# features_out = FeatureSelection(dt = dt, middle_dt = NULL,
#                                 factor = factor)
# dt = features_out$dt_new
# classifier = c("LDA", "Logit", "RF", "NNLS")
# keep_classifier = F
# features_selected = features_out$features_selected
# select_n = features_out$select_n
# test_ind = NULL
# 
# test_out = SuperSigClassifier(dt = dt, test_ind = NULL,
#                               factor = factor,
#                               classifier = c("LDA", "Logit", "RF", "NNLS"),
#                               keep_classifier = F,
#                               features_selected = features_selected,
#                               select_n = select_n)

# Create test objects for Shiny app
# test_model = SuperSigClassifier(dt = dt, test_ind = NULL,
#                                 factor = factor,
#                                 classifier = c("LDA"),
#                                 keep_classifier = T,
#                                 features_selected = features_selected)$classifier
# test_formula = features_out$features_gmsa$new_partition_formula[features_selected]
# saveRDS(list(model = test_model, formula = test_formula), "app/data/test_shiny.rds")
# write.csv(signature_caf[["Data", ind]]$DataSetFiltered[1, ], "app/data/dt.csv", row.names = F)