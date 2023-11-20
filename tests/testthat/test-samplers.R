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
      res <-.checkNamesMatchParams(names[i], as.list(rep(0, params[i])))
      
    }
    expect_true(all(res == c(d_uv, c_mv)))
  }
  
  expect_error(.checkNamesMatchParams("asld", 2))
  expect_error(.checkNamesMatchParams("norm", rep(0, 3)))
})


test_that(".checkStart", {
  info <- c(1,2)
  expect_error(.checkStart(info = c(1,1), 1))
  expect_error(.checkStart(info = c(1,0), 2))
  expect_no_error(.checkStart(info = c(1,1), 2))
})
