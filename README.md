
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplr

<!-- badges: start -->
<!-- badges: end -->

The goal of samplr is to provide tools to understand human performance
from the perspective of sampling, both looking at how people generate
samples and how people use the samples they have generated. A longer
overview and other resources can be found at
[sampling.warwick.ac.uk](sampling.warwick.ac.uk)

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
chain <- sampler_mcmc(distr_name = "norm", distr_params = c(0,1), start = 1, sigma_prop = diag(1) * .5, iterations = 2048)
print(chain[[1]][1:20])
#>  [1]  1.00000000  0.55703026  0.68688570  0.09600704  1.22404092  0.13513180
#>  [7] -0.52146425 -0.72986307 -0.73394108  0.96640564  1.50634775  1.78200855
#> [13]  1.57732831  1.57732831  1.75567722  1.75567722  0.88060541  0.72202406
#> [19]  0.98888308  1.08316613
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
