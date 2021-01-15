# all_my_auc.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Dec 10, 2019

#' Calculate AUC for every column
#' 
#' all_my_auc calls my_auc to calculate the AUC for every column other than the 
#' indicator variable
#' 
#' @param z a matrix
#' @param IndVar indicator variable for the exposure
#'
#' @return \code{all_my_auc} returns a list of AUCs, 
#' one for every column (feature)
#' 
#' @noRd
#' 
all_my_auc <- function(z, IndVar){
  cols <- setdiff(colnames(z), IndVar)
  apply(z[, cols], MARGIN = 2, FUN = my_auc, y = z %>% pull(!!IndVar))
}