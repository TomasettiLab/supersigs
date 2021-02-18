# condense_mutations.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Dec 10, 2019
#
# Helper function for GenerateMinSigmaAlgebra.R

#' Given a set of level 3 features, simplify and condense into level 1 and 
#' level 2 mutations
#' 
#' @param features a vector of level 3 features
#'
#' @return \code{condense_mutations} returns a vector of condensed features
#' 
#' @noRd
#' 
condense_mutations = function(features){
    out <- vector('character')
    for(i in seq_along(muts_children_level3)){
        v <- muts_children_level3[[i]]
        if(all(v %in% features)){
            out <- c(out, names(muts_children_level3)[i])
            features <- setdiff(features, v)
        }
    }
    return(out)
}
