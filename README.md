
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supersigs

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- [![CRAN Status](https://www.r-pkg.org/badges/version/pkgdown)](https://cran.r-project.org/package=pkgdown) -->
<!-- [![R build status](https://github.com/r-lib/pkgdown/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/supersigs/actions) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/r-lib/pkgdown/branch/master/graph/badge.svg)](https://codecov.io/gh/r-lib/supersigs?branch=master) -->
<!-- badges: end -->

`supersigs` is a companion R package to a method proposed by *Afsari, et
al. (2021, ELife)* to generate mutational signatures from single
nucleotide variants in the cancer genome. **Note: Package is under
active development.**

More details on the statistical method can be found in this paper:

-   Afsari, B., Kuo, A., Zhang, Y., Li, L., Lahouel, K., Danilova, L.,
    Favorov, A., Rosenquist, T. A., Grollman, A. P., Kinzler, K. W.,
    Cope, L., Vogelstein, B., & Tomasetti, C. (2021). Supervised
    mutational signatures for obesity and other tissue-specific
    etiological factors in cancer. ELife, 10.
    [https://doi.org/10.7554/elife.61082](https://doi.org/10.7554/eLife.61082)

## Installation

You can install the development version of supersigs from github using
the `install_github()` function from the `devtools` package.

``` r
# Install development version from GitHub
devtools::install_github("TomasettiLab/supersigs")
```

## Data format

At a minimum, the data you will need are the age and mutations for each
patient. An example is provided below. (Note that you will need to
process the data before running the core functions, see
`vignette("supersigs")` for details.)

    #>   sample_id age chromosome  position ref alt
    #> 1         1  50       chr1  94447621   G   C
    #> 2         1  50       chr2 202005395   A   C
    #> 3         1  50       chr7  20784978   T   A
    #> 4         1  50       chr7  87179255   C   G
    #> 5         1  50      chr19   1059712   G   T
    #> 6         2  55       chr1  76226977   T   C

## Core functions

In brief, the `supersigs` package contains three core functions:
`get_signature`, `predict_signature`, and `partial_signature`.

`get_signature` trains a supervised signature for a given factor
(e.g. smoking).

``` r
supersig <- get_signature(data = data, factor = "smoking", wgs = F)
```

`predict_signature` uses the trained supervised signature to obtain
predicted probabilities (e.g. probability of smoker) on a new dataset.

``` r
pred <- predict_signature(object = supersig, newdata = data, factor = "smoking")
```

`partial_signature` removes the contribution of a trained signature from
the dataset.

``` r
data <- partial_signature(data = data, object = supersig)
```

## Tutorial

To follow a tutorial on how to use the package, see
`vignette("supersigs")` (or type `vignette("supersigs")` in R).
