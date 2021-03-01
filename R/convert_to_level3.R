# convert_to_level3.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Mar 1, 2021
#
# Helper function for GenerateMinSigmaAlgebra.R

#' Given a feature, convert it to a set of 
#' fundamental (level 3) mutation features
#' 
#' @param feature a feature (e.g. "C>T")
#' 
#' @return \code{convert_to_level3} returns a vector of 
#' fundamental (level 3) features
#' 
#' @noRd
#' 
convert_to_level3 = function(feature){
    out <- unname(muts_children_level3[[feature]])
    return(out)
}
