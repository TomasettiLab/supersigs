# globals.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 11, 2021
#
# Fix visible binding note

#' @import utils
if(getRversion() >= "2.15.1"){  
  utils::globalVariables(c(".", ":=", "background_probs_wgs"))
}