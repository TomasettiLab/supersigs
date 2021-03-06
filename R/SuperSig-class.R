# SuperSig-class.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Feb 25, 2021
#' An S4 class for SuperSig
#'
#' @slot Signature data frame of features and their 
#' difference in mean rates between exposed and unexposed
#' (or the average rate if the factor is "age")
#' @slot Features list of features that comprise the signature and 
#' their representation in terms of the fundamental (trinucleotide) mutations
#' @slot AUC length-one numeric vector of the apparent AUC 
#' (i.e. not cross-validated)
#' @slot Model list of a glm class for trained logistic regression model
#' 
SuperSig <- setClass("SuperSig",
                     slots = list(Signature = "data.frame",
                                  Features = "list",
                                  AUC = "numeric",
                                  Model = "list"),
                     contains = "glm")


# Display ---
# class method for show() generic
setMethod("show",
                    "SuperSig",
                    function(object) {
                        cat("Signature:\n")
                        print(object@Signature)
                        cat("Features:\n")
                        print(object@Features)
                        cat("Model:\n")
                        print(object@Model)
                    }
)

# class method for print
print.SuperSig = function(object){
    object
}

# Accessors ---
setGeneric("Signature", function(x) standardGeneric("Signature"))
setMethod("Signature", "SuperSig", function(x) x@Signature)

setGeneric("Features", function(x) standardGeneric("Features"))
setMethod("Features", "SuperSig", function(x) x@Features)

setGeneric("AUC", function(x) standardGeneric("AUC"))
setMethod("AUC", "SuperSig", function(x) x@AUC)

setGeneric("Model", function(x) standardGeneric("Model"))
setMethod("Model", "SuperSig", function(x) x@Model)