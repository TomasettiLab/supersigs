# SaveExternalData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 7, 2020
# Save rds files to external data for R package

# Generate example VCF data for vignette
example_dt = data.frame(sample_id = rep(1:5, each = 5),
                        age = rep(c(50, 55, 72, 53, 48), each = 5),
                        chromosome = c("chr1", "chr2", "chr7", "chr7", "chr19",
                                       "chr1", "chr3", "chr9", "chr16", "chr22",
                                       "chr1", "chr4", "chr2", "chr11", "chr2",
                                       "chr1", "chr18", "chr6", "chr6", "chr1",
                                       "chr1", "chr4", "chr19", "chr6", "chr10"),
                        position = c(94447621, 202005395, 20784978, 87179255, 1059712,
                                     76226977, 38180872, 139905080, 1562631, 42189307,
                                     94447621, 202005395, 20784978, 87179255, 1059712,
                                     76226977, 38180872, 139905080, 1562631, 42189307,
                                     94447621, 202005395, 20784978, 87179255, 1059712,),
                        from = c("G", "A", "T", "C", "G",
                                 "T", "C", "G", "G", "C",
                                 "T", "A", "T", "A", "T",
                                 "G", "C", "G", "G", "G",
                                 "A", "C", "T", "G", "G"),
                        to = c("C", "C", "A", "G", "T",
                               "C", "G", "T", "T", "A",
                               "C", "C", "A", "G", "T",
                               "C", "G", "T", "T", "A",
                               "C", "C", "A", "G", "T"))

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

# Save external data
usethis::use_data(example_dt, transform_muts_vec)
