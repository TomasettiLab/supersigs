# MyAuc.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Jul 15, 2019
#
# Functions for calculating AUC using the Mann-Whitney U Statistic

# source(here("code", "AllMyAuc.R"))

#' Calculate AUC for one column using the Mann-Whitney U statistic
#' 
#' @param y Indicator variable column, i.e. "true value"
#' @param x A feature column
#' 
#' @return output AUC value
#' 
MyAuc <- function(y, x){
  # Handle missing observations
  missing_x = is.na(x)
  if(all(missing_x)){
    return(NA)
  } else if(any(missing_x)){
    x = x[!missing_x]
    y = y[!missing_x]
  }
  ind <- which(y == 1)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  diff_rank <- sum(rank(x)[ind]) - n1*(n1+1)/2 # diff_rank = 0 if all the largest values in x were y == 1
  return(diff_rank/(n1*n0))                    # n1*n0 = max difference possible
}
