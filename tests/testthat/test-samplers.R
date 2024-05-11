test_that(".checkMVInputs", {
  v <- c(0, 1)
  m <- diag(2)
  n <- 1
  
  expect_no_error(.checkMVInputs("mvnorm", list(v,m)))
  expect_error(.checkMVInputs("mvnorm", list(v,v)))
  expect_error(.checkMVInputs("mvnorm", list(m,m)))
  
  expect_no_error(.checkMVInputs("mvt", list(v,m,n)))
  expect_error(.checkMVInputs("mvt", list(m,m,n)))
  expect_error(.checkMVInputs("mvt", list(v,v,n)))
  expect_error(.checkMVInputs("mvt", list(v,m,m)))
})
  
test_that(".checkNamesMatchParams", {
  names_cont <- c("unif", "norm","lnorm", "gamma", "beta", "nbeta", "chisq", "nchisq", "t", "nt", "f", "nf", "cauchy", "exp", "logis", "weibull",
                  "4beta", "lst", "truncnorm", "trunct", "trunclst", "triangular")
  names_cont_mv <- c("mvnorm", "mvt")
  names_discr <- c("binom", "nbinom", "nbinom_mu", "pois", "geom", "hyper", "wilcox", "signrank")
  names <- c(names_cont, names_discr, names_cont_mv)
  parameters_cont <- c(2, 2, 2, 2, 2, 3, 1, 2, 1, 2, 2, 3, 2, 1, 2, 2,
                       4, 3, 4, 3, 5, 3)
  parameters_discr <- c(2,2,2,1,1,3,2,1)
  parameters_cont_mv <- c(2, 3)
  params <- c(parameters_cont, parameters_discr, parameters_cont_mv)
  
  
  # data.frame(names, params)
  for (i in 1:length(names)){
    c_uv = is.element(names[i], names_cont)
    c_mv = is.element(names[i], names_cont_mv);
    d_uv = is.element(names[i], names_discr);
    
    if (!(names[i] %in% names_cont_mv)){
      res <-.checkNamesMatchParams(names[i], rep(0, params[i]))
    } else{
      if (names[i] == "mvnorm"){
        res <-.checkNamesMatchParams(
          names[i], list(c(0,0), diag(2)))  
      } else {
        if (names[i] == "mvt"){
          res <-.checkNamesMatchParams(
            names[i], list(c(0,0), diag(2), 1))  
        }
      }
    }
    expect_true(all(res == c(d_uv, c_mv)))
  }
  
  expect_error(.checkNamesMatchParams("asld", 2))
  expect_error(.checkNamesMatchParams("norm", rep(0, 3)))
  expect_error(.checkNamesMatchParams("norm", list(rep(0, 2))))
})

test_that(".checkStart", {
  info <- c(1,2)
  expect_error(.checkStart(info = c(1,1), 1))
  expect_error(.checkStart(info = c(1,0), 2))
  expect_no_error(.checkStart(info = c(1,1), 2))
})

test_that(".checkWeights", {
  expect_warning(.checkWeights(1,1))
  expect_error(.checkWeights(1,2))
  expect_error(.checkWeights(c(1,1), 2))
  expect_message(.checkWeights(NULL, 3))
})

test_that(".checkSigmaProp", {
  expect_warning(.checkSigmaProp(NULL, 1))
  expect_warning(.checkSigmaProp(NULL, 2))
  expect_true(is.matrix(.checkSigmaProp(c(1), 1)))
  expect_no_error(.checkSigmaProp(diag(2), 2))
})

test_that(".checkGivenInfo", {
  expect_warning(.checkGivenInfo(
    distr_name="norm", 
    distr_params=c(0,1), 
    start=1, 
    weights=NULL, 
    caller="mh", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_no_warning(.checkGivenInfo(
    distr_name="norm", 
    distr_params=c(0,1), 
    start=1, 
    weights=NULL, 
    caller="hmc", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_warning(.checkGivenInfo(
    distr_name="norm", 
    distr_params=c(0,1), 
    start=1, 
    weights=NULL, 
    caller="hmc", 
    custom_density=\(x){F}, 
    sigma_prop = NULL
  ))
  expect_error(.checkGivenInfo(
    distr_name=NULL,
    distr_params=NULL, 
    start=1, 
    weights=1, 
    caller="hmc", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_error(.checkGivenInfo(
    distr_name=c("norm", "norm"),
    distr_params=list(c(1,2)), 
    start=1, 
    weights=1, 
    caller="hmc", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_error(.checkGivenInfo(
    distr_name=c("norm", "binom"),
    distr_params=list(c(1,2),c(1,2)), 
    start=1, 
    weights=1, 
    caller="hmc", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_no_error(.checkGivenInfo(
    distr_name=c("norm", "norm"),
    distr_params=list(c(1,2),c(0,2)), 
    start=1, 
    weights=c(.5, .5), 
    caller="hmc", 
    custom_density=NULL, 
    sigma_prop = NULL
  ))
  expect_warning(.checkGivenInfo(
    distr_name=NULL,
    distr_params=NULL,
    start=1, 
    weights=c(.5, .5), 
    caller="hmc", 
    custom_density=\(x){dnorm(x)}, 
    sigma_prop = NULL
  ))
})

test_that("mh", {
  expect_no_error(
    sampler_mh(0, "norm", c(0,1), diag(1))
  )
  expect_no_error(
    sampler_mh(0, sigma_prop =  diag(1), custom_density = \(x){dnorm(x)})
  )
})

test_that("mc3", {
  expect_no_error(
    sampler_mc3(start = matrix(c(0,0),ncol=1), distr_name = "norm", distr_params = c(0,1), sigma_prop = diag(1), nChains = 2)
  )
  expect_no_error(
    sampler_mc3(start = matrix(c(0,0),ncol=1), custom_density = \(x){dnorm(x)}, sigma_prop = diag(1), nChains = 2)
  )
  expect_no_error(
    sampler_mc3(start = 0, custom_density = \(x){dnorm(x)}, sigma_prop = diag(1), nChains = 2)
  )
  expect_error(
    sampler_mc3(start = 0, distr_name = "norm", distr_params = c(0,1), sigma_prop = diag(1), nChains = 2.5)
  )
  expect_error(
    sampler_mc3(start = matrix(c(0,0,0),ncol=1), distr_name = "norm", distr_params = c(0,1), sigma_prop = diag(1), nChains = 2)
  )
})

test_that("hmc", {
  expect_no_error(
    sampler_hmc(0, "norm", c(0,1))
  )
  expect_no_error(
    sampler_hmc(0, custom_density = \(x){dnorm(x)})
  )
  
  expect_error(
    sampler_hmc(0, "binom", c(0,1))
  )
})

test_that("rec", {
  expect_no_error(
    sampler_rec(0, "norm", c(0,1))
  )
  expect_no_error(
    sampler_rec(0, custom_density = \(x){dnorm(x)})
  )
  expect_error(
    sampler_rec(0, "binom", c(0,1))
  )
})

test_that("mchmc", {
  expect_no_error(
    sampler_mchmc(0, "norm", c(0,1))
  )
  expect_no_error(
    sampler_mchmc(0, custom_density = \(x){dnorm(x)})
  )
  expect_error(
    sampler_mchmc(0, "norm", c(0,1), nChains = 2.5)
  )
  expect_error(
    sampler_mchmc(0, "binom", c(0,1))
  )
})

test_that("sampler_mcrec", {
  expect_no_error(
    sampler_mcrec(0, "norm", c(0,1))
  )
  expect_no_error(
    sampler_mcrec(0, custom_density = \(x){dnorm(x)})
  )
  expect_error(
    sampler_mcrec(0, "norm", c(0,1), nChains = 2.5)
  )
  expect_error(
    sampler_mcrec(0, "binom", c(0,1))
  )
})
