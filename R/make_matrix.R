# make_matrix.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Dec 23, 2020
#
# (Export) Function to transform data frame of mutations into correct format

# library(dplyr)
# library(tidyr)

#' Function to transform mutations into "matrix" format
#' 
#' Transform a data frame of mutations in VCF format into 
#' a data frame of trinucleotide mutations with flanking bases
#' in a wide matrix format.
#' 
#' @param dt a data frame of mutations in VCF format (see vignette for details)
#' @param genome the reference genome used ("hg19" or "hg38")
#' 
#' @import dplyr
#' @import tidyr
#' 
#' @export
#' 
#' @return \code{make_matrix} returns a data frame of mutations
#' 
#' @examples
#' 
#' # print(example_dt) # use example data from package
#' # input_dt = make_matrix(example_dt) # convert to correct format
#' # get_signature(dt = input_dt, factor = "age") # get SuperSig
#' 
make_matrix <- function(dt, genome = "hg19"){
  dt = dt %>%
    select(sample_id, age, chromosome, position, from, to) %>%
    mutate(start = position - 1,
           end = position + 1)
  
  if(genome == "hg19"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      stop("Package \"BSgenome.Hsapiens.UCSC.hg19\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    dt_ranges <- as(dt %>% select(chromosome, start, end), "GRanges")
    aligned_dna <- getSeq(BSgenome.Hsapiens.UCSC.hg19, dt_ranges)
  } else if(genome == "hg38"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("Package \"BSgenome.Hsapiens.UCSC.hg38\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    dt_ranges <- as(dt %>% select(chromosome, start, end), "GRanges")
    aligned_dna <- getSeq(BSgenome.Hsapiens.UCSC.hg38, dt_ranges)
  } else {
    stop("Invalid genome specified")
  }
  
  # Create mutations with surrounding base pairs
  dt <- dt %>%
    mutate(aligned = as.character(aligned_dna),
           mutation = paste0(substr(aligned, 1, 1), "[", from, ">", to, "]",
                             substr(aligned, 3, 3)),
           mutation_std = unname(sapply(mutation, function(x) transform_muts_vec[[x]])))
  
  # Count mutations for each patient
  dt_counts <- dt %>%
    group_by(sample_id, age, mutation_std) %>%
    summarize(mut_count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = mutation_std, values_from = mut_count) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  # Add any fundamental mutations that are missing
  for(mut in transform_muts_vec){
    if(!(mut %in% names(dt_counts))){
      dt_counts = dt_counts %>%
        mutate(!!mut := 0)
    }
  }
  
  return(dt_counts)
}