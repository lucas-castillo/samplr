## R CMD check results

1 error | 0 warnings | 6 notes


❯ checking tests ...
  See below...

❯ checking C++ specification ... NOTE
    Specified C++11: please drop specification unless essential

❯ checking package subdirectories ... NOTE
  Problems with news in 'NEWS.md':
  No news entries found.

❯ checking dependencies in R code ... NOTE
  Namespaces in Imports field not imported from:
    'R6' 'gridExtra' 'magrittr' 'tibble' 'tidyr'
    All declared Imports should be used.

❯ checking R code for possible problems ... NOTE
  Bayesian_Sampler: no visible global function definition for 'sd'
  Mean_Variance: no visible binding for global variable 'var'
  Mean_Variance: no visible global function definition for 'coef'
  calc_PSD: no visible global function definition for 'abline'
  calc_all: no visible global function definition for 'par'
  calc_autocorr: no visible global function definition for 'lines'
  calc_levy: no visible global function definition for 'abline'
  calc_levy: no visible global function definition for 'lm'
  calc_qqplot: no visible global function definition for 'qqnorm'
  calc_qqplot: no visible global function definition for 'qnorm'
  calc_qqplot: no visible global function definition for 'quantile'
  calc_qqplot: no visible global function definition for 'abline'
  calc_sigma_scaling: no visible global function definition for 'abline'
  plot_2d_density: no visible binding for global variable 'x'
  plot_2d_density: no visible binding for global variable 'y'
  Undefined global functions or variables:
    abline coef lines lm par qnorm qqnorm quantile sd var x y
  Consider adding
    importFrom("graphics", "abline", "lines", "par")
    importFrom("stats", "coef", "lm", "qnorm", "qqnorm", "quantile", "sd",
               "var")
  to your NAMESPACE file.

❯ checking Rd line widths ... NOTE
  Rd file 'Zhu23ABS.Rd':
    \examples lines wider than 100 characters:
       zhuabs <- Zhu23ABS$new(nd_time = 0.3, s_nd_time = 0.5, lambda = 10, n_chains = 5, width = 1, distr_name = 'norm')
       tafc_sim <- zhuabs$two_alt_force_choice(delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim)
  
  Rd file 'sampler_mchmc.Rd':
    \examples lines wider than 100 characters:
       MCHMC <- sampler_mchmc(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)
  
  Rd file 'sampler_mcrec.Rd':
    \examples lines wider than 100 characters:
       MCREC <- sampler_mcrec(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)
  
  These lines will be truncated in the PDF manual.

❯ checking files in 'vignettes' ... NOTE
  Files named as vignettes but with no recognized vignette engine:
     'vignettes/ABS-Parameter-Recovery.Rmd'
  (Is a VignetteBuilder field missing?)
```
── Test failures ─────────────────────────────────────────────────────────────────────────────────────── testthat ────

> # This file is part of the standard setup for testthat.
> # It is recommended that you do not modify it.
> #
> # Where should you do additional test configuration?
> # Learn more about the roles of various files in:
> # * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
> # * https://testthat.r-lib.org/articles/special-files.html
> 
> library(testthat)
> library(samplr)
> 
> test_check("samplr")
[ FAIL 5 | WARN 0 | SKIP 0 | PASS 173 ]

══ Failed tests ════════════════════════════════════════════════════════════════
── Failure ('test-calc_functions.R:77:3'): Sigma Scaling ───────────────────────
`res` has type 'list', not 'double'.
── Failure ('test-calc_functions.R:81:3'): Sigma Scaling ───────────────────────
`res` (`actual`) not equal to c(16.40681, 20.13771, 19.16075, 17.64051, 15.55798, 17.40862) (`expected`).

`actual` is a list
`expected` is a double vector (16.40681, 20.13771, 19.16075, 17.64051, 15.55798, ...)
── Error ('test-calc_functions.R:91:3'): Autocorr Plotter ──────────────────────
Error in `plot_autocorr(matrix(sequence, ncol = 3))`: could not find function "plot_autocorr"
Backtrace:
    ▆
 1. └─testthat::expect_error(...) at test-calc_functions.R:91:3
 2.   └─testthat:::expect_condition_matching(...)
 3.     └─testthat:::quasi_capture(...)
 4.       ├─testthat (local) .capture(...)
 5.       │ └─base::withCallingHandlers(...)
 6.       └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
── Failure ('test-calc_functions.R:96:3'): Series Plotter ──────────────────────
`plot_series(matrix(sequence, ncol = 3))` did not throw the expected error.
── Error ('test-calc_functions.R:101:3'): Change Plotter ───────────────────────
Error in `plot_change(matrix(sequence, ncol = 3))`: could not find function "plot_change"
Backtrace:
    ▆
 1. └─testthat::expect_error(...) at test-calc_functions.R:101:3
 2.   └─testthat:::expect_condition_matching(...)
 3.     └─testthat:::quasi_capture(...)
 4.       ├─testthat (local) .capture(...)
 5.       │ └─base::withCallingHandlers(...)
 6.       └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))

[ FAIL 5 | WARN 0 | SKIP 0 | PASS 173 ]
Error: Test failures
Execution halted
```