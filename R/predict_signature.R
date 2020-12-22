# predict_signature.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 22, 2020
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
#' 
#' @export
#' 
#' @return \code{predict_signature} returns the original data frame with 
#' additional columns for the feature counts and classification score
#' 
#' @examples
#' 
#' 
predict_signature <- function(object, newdata, factor){
  # Extract slots from object
  model = object@Model$Logit
  features = object@Features
  
  # Add counts for features that are used in the model
  for(j in seq_along(features)){
    feature = names(features)[j]
    feature_count = newdata %>%
      select(features[[feature]]) %>%
      rowSums()
    newdata = newdata %>%
      mutate(!!feature := feature_count)
  }
  
  # Use rates for non-age factors
  if(factor != "AGE"){
    newdata = newdata %>%
      mutate_at(names(features), ~(./age))
  }
  
  # Predict using logistic regression
  scores = predict(model, 
                   newdata = newdata, 
                   type = "response")
  
  newdata = newdata %>%
    mutate(score = scores)
  
  return(newdata)
}