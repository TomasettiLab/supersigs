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
#' in a wide matrix format.
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
  sample_ids <- rownames(colData(vcf))
  if("age" %in% colnames(colData(vcf))){
    ages <- colData(vcf)$age
  } else {
    ages <- rep(NA, n = n_samples)
  }
  
  # positions for homozygous or heterozygous alt
  positions <- geno(vcf)$GT != "0|0"
  
  dt_ls <- vector("list", n_samples)
  for(i in seq_len(n_samples)){
    dt_ls[[i]] <- data.frame(rowRanges(vcf[positions[, i],])) %>%
      dplyr::mutate(chromosome = paste0("chr", seqnames)) %>%
      dplyr::rename(position = start) %>%
      dplyr::mutate(ref = vapply(REF, function(x)
      {as.character(x)[[1]][[1]]},
      FUN.VALUE = character(1)),
      alt = vapply(ALT, function(x)
      {as.character(x)[[1]][[1]]},
      FUN.VALUE = character(1)),
      sample_id = sample_ids[i],
      age = ages[i]) %>%
      dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1) %>%
      dplyr::select(sample_id, age, chromosome, position, ref, alt)
  }
  dt <- rbind(dt_ls[[i]])
  
  return(dt)
}
