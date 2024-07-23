
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplr: Tools To Compare Human Performance To Sampling Algorithms

<!-- badges: start -->


[![R-CMD-check](https://github.com/lucas-castillo/samplr/workflows/R-CMD-check/badge.svg)](https://github.com/lucas-castillo/samplr/actions)
[![Codecov test
coverage](https://codecov.io/gh/lucas-castillo/samplr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lucas-castillo/samplr?branch=main)
[![License: CC BY
4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R-CMD-check](https://github.com/lucas-castillo/samplr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lucas-castillo/samplr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of samplr is to provide tools to understand human performance
from the perspective of sampling, both looking at how people generate
samples and how people use the samples they have generated. A longer
overview and other resources can be found at
[sampling.warwick.ac.uk](sampling.warwick.ac.uk). Get started
[here](vignettes/how-to-sample.html).

## Installation

You can install samplr from CRAN:

    install.packages("samplr")

or install the development version from
[Github](https://github.com/lucas-castillo/samplr) with:

    devtools::install_github("lucas-castillo/samplr")

or alternatively using the `remotes` package

    remotes::install_github("lucas-castillo/samplr")

### Installing development version on MacOS

If you are installing the development version on MacOS, you will need
the following prior to installation:

1.  Apple’s ‘Command Line Tools’: these can be (re-)installed by running
    `xcode-select --install` in a terminal. You may also check if those
    are already installed by running `pkgbuild::check_build_tools()` in
    R.
2.  A Fortran compiler. Installers for gfortran are available
    [here](https://github.com/fxcoudert/gfortran-for-macOS/releases/).
    This installs into `/usr/local/gfortran`.

Read more about it on the [macOS Prerequisites section in the R
Installation and Administration
Manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Prerequisites).

### Installing development version on Windows

If you are installing the development version on Windows, you will need
to have RTools installed, which you can find
[here](https://cran.r-project.org/bin/windows/Rtools/). Please make sure
you install the version corresponding to your R version (i.e. for R
4.3.3, you’d need RTools 4.3).

## Example

The samplr package provides tools to generate samples following
particular algorithms

``` r
library(samplr)
set.seed(1)
chain <- sampler_mh(start = 1, distr_name = "norm", distr_params = c(0,1), sigma_prop = diag(1) * .5, iterations = 2048)
r <- plot_series(chain[[1]], change = FALSE)
```

![](man/figures/README-unnamed-chunk-2-1.png)<!-- -->

As well as tools to diagnose the patterns both from samplers and
participants:

``` r
v <- calc_all(chain[[1]][1:200])
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->
<!-- TODO: add plot_series(chain[[1]]) again -->
