# SaveInternalData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 11, 2021
# Save rds files to internal data for R package

library(here)
library(devtools)

# Data created by Bahman

all_possible_tri = readRDS(here("data-files", "all_possible_tri.rds"))
background_probs = readRDS(here("data-files", "background_probs.rds"))
background_probs_wgs = readRDS(here("data-files", "background_probs_wgs.rds"))
h_muts_index = readRDS(here("data-files", "h_muts_index.rds"))
muts_children_level3_df = readRDS(here("data-files", "muts_children_level3_df.rds"))
muts_children_level3 = readRDS(here("data-files", "muts_children_level3.rds"))
muts_formula = readRDS(here("data-files", "muts_formula.rds"))
muts_level3 = readRDS(here("data-files", "muts_level3.rds"))
prop_muts_all = readRDS(here("data-files", "prop_muts_all.rds"))

# Create transform_muts_vec, a named vector that converts mutations to a standard format,
# i.e. its fundamental mutation must begin with "C" or T"
StandardizeMutations <- function(){
  base_pair <- setNames(c("C", "T", "A", "G"), c("G", "A", "T", "C"))
  mutation_vec <- c()
  mutation_std_vec <- c()
  
  for(orig in base_pair){
    for(new in base_pair){
      if(orig == new) next
      for(begin in base_pair){
        for(end in base_pair){
          mutation <- paste0(base_pair[begin], "[", 
                             base_pair[orig], ">", base_pair[new], "]",
                             base_pair[end])
          mutation_complement <- paste0(names(base_pair[end]), "[", 
                                        names(base_pair[orig]), ">", 
                                        names(base_pair[new]), "]",
                                        names(base_pair[begin]))
          mutation_vec <- c(mutation_vec, mutation)
          if(!(orig %in% c("C", "T"))){
            mutation_std_vec <- c(mutation_std_vec, mutation)
          } else {
            mutation_std_vec <- c(mutation_std_vec, mutation_complement)
          }
        }
      }
    }
  }
  
  out <- setNames(mutation_std_vec, mutation_vec)
  return(out)
}

transform_muts_vec <- StandardizeMutations() # named vector

# Write internal data to R/sysdata.rda
use_data(all_possible_tri,
         background_probs,
         h_muts_index,
         muts_children_level3_df,
         muts_children_level3,
         muts_formula,
         muts_level3,
         prop_muts_all, 
         transform_muts_vec,
         internal = TRUE,
         overwrite = TRUE)