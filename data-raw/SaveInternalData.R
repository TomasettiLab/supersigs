# SaveInternalData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 10, 2019
# Save rds files to internal data for R package

library(here)
library(devtools)

all_possible_tri = readRDS(here("data", "all_possible_tri.rds"))
background_probs = readRDS(here("data", "background_probs.rds"))
h_muts_index = readRDS(here("data", "h_muts_index.rds"))
muts_children_level3_df = readRDS(here("data", "muts_children_level3_df.rds"))
muts_children_level3 = readRDS(here("data", "muts_children_level3.rds"))
muts_formula = readRDS(here("data", "muts_formula.rds"))
muts_level3 = readRDS(here("data", "muts_level3.rds"))
prop_muts_all = readRDS(here("data", "prop_muts_all.rds"))

# Write internal data to R/sysdata.rda
use_data(all_possible_tri,
         background_probs,
         h_muts_index,
         muts_children_level3_df,
         muts_children_level3,
         muts_formula,
         muts_level3,
         prop_muts_all, internal = TRUE)