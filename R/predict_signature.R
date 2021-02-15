# predict_signature.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 23, 2020
#
# (Export) Function to apply SuperSigs model on new dataset

#' Function to predict using SuperSig object
#' 
#' Using a generated SuperSig, predict on a new dataset and return 
#' predicted probabilities for each observation. 
#' 
#' @param object an object of class "SuperSig" 
#' @param newdata a data frame of mutations
#' @param factor the factor/exposure (e.g. "age", "smoking")
#' 
#' @import dplyr
#' @importFrom methods as new
#' @importFrom stats binomial coef end glm median pbinom predict setNames start
#' 
#' @export
#' 
#' @return \code{predict_signature} returns the original data frame with 
#' additional columns for the feature counts and classification score
#' 
#' @examples
#' 
#' head(example_dt) # use example data from package
#' input_dt <- make_matrix(example_dt) # convert to correct format
#' input_dt$IndVar <- c(1, 1, 1, 0, 0) # add IndVar column
#' out <- get_signature(dt = input_dt, factor = "age") # get SuperSig
#' 
#' newdata <- predict_signature(out, newdata = input_dt, factor = "Smoking")
#' suppressPackageStartupMessages({library(dplyr)})
#' head(newdata %>% select(IndVar, F1, score))
#' 
predict_signature <- function(object, newdata, factor){
  # Extract slots from object
  model <- Model(object)$Logit
  features <- Features(object)
  
  # Add counts for features that are used in the model
  for(j in seq_along(features)){
    feature <- names(features)[j]
    feature_count <- newdata %>%
      select(features[[feature]]) %>%
      rowSums()
    newdata <- newdata %>%
      mutate(!!feature := feature_count)
  }
  
  # Use rates for non-age factors
  if(factor != "AGE"){
    newdata <- newdata %>%
      mutate_at(names(features), ~(./age))
  }
  
  # Predict using logistic regression
  scores <- predict(model, 
                    newdata = newdata, 
                    type = "response")
  
  newdata <- newdata %>%
    mutate(score = scores)
  
  return(newdata)
}