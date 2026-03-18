# Development version
* Bugfixes:
   * Fixes calculation of some probabilities in Z_identities; p(A or not B) and p(B or not A)
* Minor:
    * Fixes typo in "Get started" vignette.
# samplr 1.1.1
* Bugfixes:
    * `Mean_Variance()` function now returns as many rows as IDs (used to return too many)
* Minor:
    * Adds citation information to README
    * Removes `latex2exp` dependency
    * Fixes invalid URL in samplr.Rmd
    * Substitutes http links for https

# samplr 1.1.0
* New features:
    * `Zhu23ABS()` now can return estimates when using relative stopping rule.
    * The `Bayesian_Sampler()` function can now be used to not only return expected predictions but also 
        * expected variance
        * simulated predictions.
    * `calc_PSD()` now returns a named vector in `$polyfit`, so users know what's the intercept and what's the slope
    * `calc_PSD()` now filters frequencies within a range, as done in previous literature (Gilden, 1995; Zhu et al., 2022)
* Minor:
    * Fixes documentation of some functions not including the reference list.

# samplr 1.0.1
* Fixes float comparisons using `!=` and thus not passing tests in some OS.  

# samplr 1.0.0
* Initial CRAN submission.
