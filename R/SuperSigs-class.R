# SuperSigs-class.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 22, 2020
#' An S4 class for SuperSigs
#'
#' @slot Signature data frame of features and their difference in mean rates 
#' (or counts if the factor is "age")
#' @slot Features list of features that comprise the signature and 
#' their representation in terms of the fundamental (trinucleotide) mutations
#' @slot AUC length-one numeric vector of the apparent AUC (i.e. not cross-validated)
#' @slot Model list of a glm class for the trained logistic regression model
#' 
SuperSigs <- setClass("SuperSigs",
                      slots = list(Signature = "data.frame",
                                   Features = "list",
                                   AUC = "numeric",
                                   Model = "list"),
                      contains = "glm")


# Display ---
# class method for show() generic
setMethod("show",
          "SuperSigs",
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
print.SuperSigs = function(object){
  object
}



