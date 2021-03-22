
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplr

<!-- badges: start -->
<!-- badges: end -->

The goal of samplr is to provide tools to understand human performance
from the perspective of sampling, both looking at how people generate
samples and how people use the samples they have generated. A longer
overview and other resources can be found at
[sampling.warwick.ac.uk](sampling.warwick.ac.uk). Get started
[here](vignettes/how-to-sample.html).

## Installation

You can install the released version of SampleR from
[Github](https://github.com/lucas-castillo/samplr) with:

    devtools::install_github("lucas-castillo/samplr")

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
plot_qqplot(chain[[1]], change = TRUE)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

``` r
plot_series(chain[[1]])
```

![](man/figures/README-unnamed-chunk-3-2.png)<!-- -->
