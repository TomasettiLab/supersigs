# SaveExternalData.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 11, 2020
# Save rds files to external data for R package

# Generate example VCF data for vignette
example_dt = data.frame(sample_id = rep(1:5, each = 5),
                        age = rep(c(50, 55, 72, 53, 48), each = 5),
                        chromosome = c("chr1", "chr2", "chr7", "chr7", "chr19",
                                       "chr1", "chr3", "chr9", "chr16", "chr22",
                                       "chr1", "chr2", "chr7", "chr7", "chr19",
                                       "chr1", "chr3", "chr9", "chr16", "chr22",
                                       "chr1", "chr2", "chr7", "chr7", "chr19"),
                        position = c(94447621, 202005395, 20784978, 87179255, 1059712,
                                     76226977, 38180872, 139905080, 1562631, 42189307,
                                     94447621, 202005395, 20784978, 87179255, 1059712,
                                     76226977, 38180872, 139905080, 1562631, 42189307,
                                     94447621, 202005395, 20784978, 87179255, 1059712),
                        from = c("G", "A", "T", "C", "G",
                                 "T", "C", "G", "G", "C",
                                 "T", "A", "T", "A", "T",
                                 "G", "C", "G", "G", "G",
                                 "A", "C", "T", "G", "G"),
                        to = c("C", "C", "A", "G", "T",
                               "C", "G", "T", "T", "A",
                               "C", "C", "A", "G", "A",
                               "C", "G", "T", "T", "A",
                               "C", "A", "A", "T", "T"))

# Save external data
usethis::use_data(example_dt, overwrite = TRUE)
