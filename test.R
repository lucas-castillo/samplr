X <- bench::mark(
  Samplr::sampler_mcmc("norm", c(0,1), 1, sigma_prop=1),
  Samplr::sampler_mc3("norm", c(0,1), 1, sigma_prop=1),
  Samplr::sampler_hmc("norm", c(0,1), 1),
  Samplr::sampler_nuts("norm", c(0,1), 1),
  check = FALSE,
  iterations = 100
)
print(X[,c("expression", "min", "median")])

