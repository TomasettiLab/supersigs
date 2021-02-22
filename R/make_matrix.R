# make_matrix.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
# Date last modified: Feb 22, 2021
#
# (Export) Function to transform data frame of mutations into correct format

# library(dplyr)
# library(tidyr)

# Wrapper for getSeq given specified genome (helper function)
getseq_wrapper <- function(dt, genome = "hg19"){
  if(genome == "hg19"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
      stop("Package \"BSgenome.Hsapiens.UCSC.hg19\" 
                 needed for this function to work. Please install it.",
           call. = FALSE)
    }
    dt_ranges <- as(dt %>% 
                      select(.data$chromosome, .data$start, .data$end), 
                    "GRanges")
    aligned_dna <- 
      getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
             dt_ranges)
  } else if(genome == "hg38"){
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
      stop("Package \"BSgenome.Hsapiens.UCSC.hg38\" 
                 needed for this function to work. Please install it.",
           call. = FALSE)
    }
    dt_ranges <- as(dt %>% select(.data$chromosome, .data$start, .data$end), 
                    "GRanges")
    aligned_dna <- 
      getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
             dt_ranges)
  } else {
    stop("Invalid genome specified")
  }
  return(aligned_dna)
}

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
#' @importFrom Biostrings getSeq
#' @importFrom rlang .data
#' 
#' @export
#' 
#' @return \code{make_matrix} returns a data frame of mutations,
#' one row per sample
#' 
#' @examples
#' 
#' head(example_dt) # use example data from package
#' input_dt <- make_matrix(example_dt) # convert to correct format
#' head(input_dt)
#' 
make_matrix <- function(dt, genome = "hg19"){
    dt <- dt %>%
        select(.data$sample_id, .data$age, .data$chromosome, .data$position, 
               .data$ref, .data$alt) %>%
        mutate(start = .data$position - 1,
               end = .data$position + 1)
    
    aligned_dna <- getseq_wrapper(dt, genome)
    
    # Create mutations with surrounding base pairs
    dt <- dt %>%
        mutate(aligned = as.character(aligned_dna),
               mutation = paste0(substr(.data$aligned, 1, 1), "[", 
                                 .data$ref, ">", .data$alt, "]",
                                 substr(.data$aligned, 3, 3)),
               mutation_std = vapply(.data$mutation,
                                     function(x){transform_muts_vec[[x]]},
                                     FUN.VALUE = character(1)) %>% unname())
    
    # Count mutations for each patient
    dt_counts <- dt %>%
        group_by(.data$sample_id, .data$age, .data$mutation_std) %>%
        summarize(mut_count = n()) %>%
        ungroup() %>%
        pivot_wider(names_from = .data$mutation_std, 
                    values_from = .data$mut_count) %>%
        mutate_all(~replace(., is.na(.), 0))
    
    # Add any fundamental mutations that are missing
    for(mut in transform_muts_vec){
        if(!(mut %in% names(dt_counts))){
            dt_counts <- dt_counts %>%
                mutate(!!mut := 0)
        }
    }
    
    return(dt_counts)
}
