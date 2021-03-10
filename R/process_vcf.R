# process_vcf.R
# -----------------------------------------------------------------------------
# Author: Albert Kuo
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
#' @importFrom SummarizedExperiment colData rowRanges seqinfo
#' @importFrom rlang .data
#' 
#' @export
#' 
#' @return \code{process_vcf} returns a data frame of mutations, 
#' one row per mutation
#' 
#' @examples
#' 
#' # Use example vcf from VariantAnnotation
#' suppressPackageStartupMessages({library(VariantAnnotation)})
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(fl, "hg19") 
#' 
#' # Subset to first sample
#' vcf <- vcf[, 1]
#' # Subset to row positions with homozygous or heterozygous alt
#' positions <- geno(vcf)$GT != "0|0" 
#' vcf <- vcf[positions[, 1],]
#' colData(vcf)$age <- 50        # Add patient age to colData (optional)
#' 
#' # Run function
#' dt <- process_vcf(vcf)
#' head(dt)
#' 
process_vcf <- function(vcf){
    genome <- unname(genome(SummarizedExperiment::seqinfo(vcf)))
    n_samples <- nrow(colData(vcf))
    assert_that(n_samples == 1, msg = "Number of samples not equal to 1 in VCF")
    
    sample_ids <- rownames(SummarizedExperiment::colData(vcf))
    if("age" %in% colnames(SummarizedExperiment::colData(vcf))){
        ages <- colData(vcf)$age
    } else {
        ages <- rep(NA, n = n_samples)
    }
    
    dt <- data.frame(SummarizedExperiment::rowRanges(vcf)) %>%
            dplyr::mutate(chromosome = paste0("chr", .data$seqnames)) %>%
            dplyr::rename(position = .data$start) %>%
            dplyr::mutate(ref = vapply(.data$REF, function(x)
            {as.character(x)[[1]][[1]]},
            FUN.VALUE = character(1)),
            alt = vapply(.data$ALT, function(x)
            {as.character(x)[[1]][[1]]},
            FUN.VALUE = character(1)),
            sample_id = sample_ids[1],
            age = ages[1]) %>%
            dplyr::filter(nchar(.data$ref) == 1 & nchar(.data$alt) == 1) %>%
            dplyr::select(.data$sample_id, .data$age, .data$chromosome, 
                          .data$position, .data$ref, .data$alt)
    
    return(dt)
}