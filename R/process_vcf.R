# process_vcf.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 11, 2021
#
# (Export) Function to transform VCF object into correct format

# library(dplyr)
# library(tidyr)

#' Function to transform VCF object into "matrix" format
#' 
#' Transform a VCF object into 
#' a data frame of trinucleotide mutations with flanking bases
#' in a wide matrix format. The function assumes that the VCF object
#' contains only one sample and that each row in rowRanges 
#' represents an observed mutation in the sample.
#' 
#' @param vcf a VCF object (from `VariantAnnotation` package)
#' 
#' @import dplyr
#' @import tidyr
#' @import VariantAnnotation
#' @importFrom Biostrings getSeq
#' @importFrom rlang .data
#' 
#' @export
#' 
#' @return \code{process_vcf} returns a data frame of mutations, 
#' one row per mutation
#' 
#' @examples
#' 
process_vcf <- function(vcf){
  genome <- unname(genome(seqinfo(vcf)))
  n_samples <- nrow(colData(vcf))
  assert_that(n_samples == 1, msg = "Number of samples not equal to 1 in VCF")
  
  sample_ids <- rownames(colData(vcf))
  if("age" %in% colnames(colData(vcf))){
    ages <- colData(vcf)$age
  } else {
    ages <- rep(NA, n = n_samples)
  }
  
  dt <- data.frame(rowRanges(vcf)) %>%
      dplyr::mutate(chromosome = paste0("chr", seqnames)) %>%
      dplyr::rename(position = start) %>%
      dplyr::mutate(ref = vapply(REF, function(x)
      {as.character(x)[[1]][[1]]},
      FUN.VALUE = character(1)),
      alt = vapply(ALT, function(x)
      {as.character(x)[[1]][[1]]},
      FUN.VALUE = character(1)),
      sample_id = sample_ids[1],
      age = ages[1]) %>%
      dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1) %>%
      dplyr::select(sample_id, age, chromosome, position, ref, alt)
  
  return(dt)
}
