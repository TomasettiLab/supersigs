
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supersigs

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- [![CRAN Status](https://www.r-pkg.org/badges/version/pkgdown)](https://cran.r-project.org/package=pkgdown) -->
<!-- [![R build status](https://github.com/r-lib/pkgdown/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/supersigs/actions) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/r-lib/pkgdown/branch/master/graph/badge.svg)](https://codecov.io/gh/r-lib/supersigs?branch=master) -->
<!-- badges: end -->

supersigs is a companion R package to a method proposed by *Afsari, et
al.*, which generates tissue-specific mutational signatures for the
cancer genome using supervised learning.

## Installation

You can install the development version of supersigs from github using
the `install_github()` function from the `devtools` package.

``` r
# Install development version from GitHub
devtools::install_github("albertkuo/supersigs")
```

## Core functions

The `supersigs` package contains two core functions: `get_signature` and
`predict_signature`.

`get_signature` trains a supervised signature for a given factor.

``` r
sig <- get_signature(dt = data, factor = "smoking", wgs = F)
```

`predict_signature` uses the trained supervised signature to obtain
predicted probabilities (e.g. probability of smoker) on a new dataset.

``` r
pred <- predict_signature(object = sig, newdata = data, factor = "smoking")
```

For more details, see `vignette("supersigs")`.
