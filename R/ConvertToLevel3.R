# ConvertToLevel3.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 10, 2019
#
# Helper function for GenerateMinSigmaAlgebra.R

#' Given a set of features, convert them to a union set of fundamental (level 3) 
#' mutation features
#' 
#' @param features a vector of features
#' 
#' @return \code{ConvertToLevel3} returns a vector of fundamental (level 3) features
#' 
#' @noRd
#' 
ConvertToLevel3 = function(features){
  if(length(features) > 1){
    out = sapply(features, function(x) muts_children_level3[[x]])
    out = unname(unlist(out))
  } else {
    out = unname(muts_children_level3[[features]])
  }
  return(out)
}