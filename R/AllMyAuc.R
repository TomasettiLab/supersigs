# AllMyAuc.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019

#' Calculate AUC for every column
#' 
#' AllMyAuc calls MyAuc to calculate the AUC for every column other than the indicator variable
#' 
#' @param z a matrix
#' @param IndVar indicator variable for the exposure
#'
#' @return \code{AllMyAuc} returns a list of AUCs, one for every column (feature)
#' 
AllMyAuc <- function(z, IndVar){
  cols = setdiff(colnames(z), IndVar)
  apply(z[, cols], MARGIN = 2, FUN = MyAuc, y = z %>% pull(!!IndVar))
}