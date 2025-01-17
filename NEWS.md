# samplr dev
* ABS now can return estimates when using relative stopping rule.
* The Bayesian Sampler function can now be used to not only return expected predictions but also expected variance
* The Bayesian Sampler function can be now be used to return simulated predictions.
* `calc_PSD` now returns a named vector in `$polyfit`, so users know what's the intercept and what's the slope

# samplr 1.0.1
* Fixes float comparisons using `!=` and thus not passing tests in some OS.  

# samplr 1.0.0
* Initial CRAN submission.
