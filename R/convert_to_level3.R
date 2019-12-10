# GenerateMinSigmaAlgebra.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 10, 2019
#
# Helper function for GenerateMinSigmaAlgebra.R

# Given a set of features, convert them to a union set of fundamental mutation features
convert_to_level3 = function(features){
  if(length(features) > 1){
    out = sapply(features, function(x) muts_children_level3[[x]])
    out = unname(unlist(out))
  } else {
    out = unname(muts_children_level3[[features]])
  }
  return(out)
}