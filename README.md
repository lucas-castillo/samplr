
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplr: Tools To Compare Human Performance To Sampling Algorithms

<!-- badges: start -->

[![R-CMD-check](https://github.com/lucas-castillo/samplr/workflows/R-CMD-check/badge.svg)](https://github.com/lucas-castillo/samplr/actions)
[![Codecov test
coverage](https://codecov.io/gh/lucas-castillo/samplr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lucas-castillo/samplr?branch=main)
[![License: CC BY
4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
<!-- badges: end -->

The goal of samplr is to provide tools to understand human performance
from the perspective of sampling, both looking at how people generate
samples and how people use the samples they have generated. A longer
overview and other resources can be found at
[sampling.warwick.ac.uk](sampling.warwick.ac.uk). Get started
[here](vignettes/how-to-sample.html).

## Installation

You can install the released version of samplr from
[Github](https://github.com/lucas-castillo/samplr) with:

    devtools::install_github("lucas-castillo/samplr")

or alternatively use the `remotes` package

    remotes::install_github("lucas-castillo/samplr")

### Installing on MacOS

If installing on MacOS, you will need the following prior to
installation:

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

## Example

samplr provides tools to generate samples following particular
algorithms

``` r
library(samplr)
set.seed(1)
chain <- sampler_mh(start = 1, distr_name = "norm", distr_params = c(0,1), sigma_prop = diag(1) * .5, iterations = 2048)
print(chain[[1]][1:20])
#>  [1]  1.00000000  0.55703026  0.68688570  0.68688570  0.91988288  0.26328684
#>  [7]  0.05488801  0.05081000  0.05081000  0.05081000 -0.76070605 -0.76070605
#> [13] -1.05168815 -1.06313640 -0.75506178 -0.75506178 -0.10524665  0.44780723
#> [19]  1.01645509  1.45473808
```

As well as tools to diagnose the patterns both from samplers and
participants:

``` r
v <- calc_qqplot(chain[[1]], change = TRUE, plot=T)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->
<!-- TODO: add plot_series(chain[[1]]) again -->
