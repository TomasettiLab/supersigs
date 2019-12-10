# AllMyAuc.R
# -----------------------------------------------------------------------------
# Author:             Bahman Afsari, Albert Kuo
# Date last modified: Jul 15, 2019
#
# Functions for calculating AUC using the Mann-Whitney U Statistic

# Calculate AUC for every column
AllMyAuc <- function(z, DpnVar){
  cols = setdiff(colnames(z), DpnVar)
  apply(z[, cols], MARGIN = 2, FUN = MyAuc, y = z %>% pull(!!DpnVar))
}