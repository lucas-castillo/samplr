test_that("bad inputs in initialising",{
  expect_error(Zhu23ABS$new(width = "1", n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1))
  expect_error(Zhu23ABS$new(width = 1, n_chains = c(1, 3), nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1))
  expect_error(Zhu23ABS$new(width = 1, n_chains = 5, nd_time = '0.3', s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1))
  expect_error(Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = c(0.5, 0.7), lambda = 10, distr_name = 'norm', distr_params = 1))
  expect_error(Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, distr_name = 'norm', distr_params = 1))
  expect_error(Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5))
  
  custom_post <- function(x){
    if (x >= -3 & x < -1){
      return(0.25 * x + 0.75)
    } else if (x >= -1 & x < 0) {
      return(-0.25 * x + 0.25)
    } else {
      return (0)
    }
  }
  expect_message(Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 5, distr_name = 'norm',
                              distr_params = 1, custom_distr = custom_post, custom_start = 0))
  
  
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = c(10, 11))
  trial_stim <- round(runif(5, 10, 50))
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim), 'The length of "distr_params" should equal to either 1 or the length of "trial_stim".')
})


test_that("run simulation twice", {
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1)
  trial_stim <- round(runif(5, 10, 50))
  zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim)
  trial_stim <- round(runif(5, 5, 10))
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim))
})

test_that("bad inputs in simulations with the fixed stopping rule", {
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1)
  trial_stim <- round(runif(5, 10, 50))
  expect_error(zhuabs$simulate(stopping_rule='fixedd', n_sample = 5, trial_stim = trial_stim), 'The stopping rule "fixedd" is not supported by ABS.')
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = '5', trial_stim = trial_stim))
  
  trial_stim[0] <- "10"
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim))
  
  trial_stim <- round(runif(5, 10, 50)) # reset trial_stim
  start_point <- runif(4, 10, 50)
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim, start_point = start_point), 'The length of "start_point" should equal to the length of "trial_stim".')
  start_point[5] <- NA
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim, start_point = start_point), 'Argument "start_point" contains NA.')
  start_point <- rep('1', 5)
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim, start_point = start_point), 'Argument "start_point" should be a numeric vector.')
  
})


test_that("bad inputs in simulations with the relative stopping rule",{
  # the relative stopping rule
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1)
  trial_stim <- factor(c('left', 'left', 'right', 'right', 'right'))
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = "4", dec_bdry = 0, discrim = 1, trial_stim = trial_stim), 'Argument "delta" should be a single integer.')
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 3, dec_bdry = 0, discrim = 1, trial_stim = trial_stim, prior_on_resp = c(3, 1)), "The relative difference in the prior on responses should be smaller than the relative stopping rule before the sampling process. Please adjust \"delta\" or \"prior_on_resp\".")
  expect_no_error(zhuabs$simulate(stopping_rule='relative', delta = 3, dec_bdry = 0, discrim = 1, trial_stim = trial_stim, prior_on_resp = c(3, 1), prior_depend = FALSE))
  zhuabs$reset_sim_results()
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = "0", discrim = 1, trial_stim = trial_stim))
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = "1", trial_stim = trial_stim))
  
  trial_stim <- factor(c('left', 'left', 'right', 'up', 'right'))
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim), "Argument \"trial_stim\" should not have more than two levels.")
  
  trial_stim <- c('left', 'left', 'right', 'right', 'right')
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim), 'Argument "trial_stim" should be a factor.')
  
  start_point <- runif(4, -3, 3)
  start_point[5] <- NA
  trial_stim <- factor(c('left', 'left', 'right', 'right', 'right'))
  expect_error(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim, start_point = start_point))
})


test_that("bad inputs in the confidence interval function",{
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, lambda = 10, distr_name = 'norm', distr_params = 1)
  trial_stim <- round(runif(5, 10, 50))
  expect_error(zhuabs$confidence_interval(0.5), "Please run the `estimate` method first.\n")
  expect_error(zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim)$confidence_interval(1.1), 'Argument "conf_level" should be a single value between 0 and 1.')
  zhuabs$reset_sim_results()
  trial_stim <- factor(c('left', 'left', 'right', 'left', 'right'))
  expect_warning(zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim)$confidence_interval(0.5))
})


test_that("starting points",{
  # the fixed stopping rule
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0, lambda = 10, distr_name = 'norm', distr_params = rep(1, 5))
  trial_stim <- round(runif(5, 10, 50))
  start_point <- runif(5, 10, 50)
  zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim, start_point = start_point)
  first_sample <- sapply(zhuabs$sim_results$samples, function(samples) samples[1])
  expect_equal(first_sample, start_point)
  
  # the relative stopping rule
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0, lambda = 10, distr_name = 'norm', distr_params = rep(1, 5))
  trial_stim <- factor(sample(c('left', 'right'), 5, TRUE))
  start_point <- runif(5, -3, 3)
  zhuabs$simulate(stopping_rule='relative', delta = 4, dec_bdry = 0, discrim = 3, trial_stim = trial_stim, start_point = start_point)
  first_sample <- sapply(zhuabs$sim_results$samples, function(samples) samples[1])
  expect_equal(first_sample, start_point)
})


test_that("the fixed stopping rule", {
  zhuabs <- Zhu23ABS$new(width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0, lambda = 10, distr_name = 'norm', distr_params = 1)
  trial_stim <- round(runif(5, 10, 50))
  n_sample = round(runif(1, 2, 20))
  zhuabs$simulate(stopping_rule='fixed', n_sample = n_sample, trial_stim = trial_stim)
  counts <- sapply(zhuabs$sim_results$samples, function(samples) length(samples))
  expect_equal(rep(n_sample, length(trial_stim)), counts)
})


test_that("The custom starting point with fixed stopping rule",{
  trial_stim <- sample(20:25, 2, replace=TRUE)
  custom_post_func <- function(x){
    if (x >= 19 & x < 22){
      return(0.3)
    } else if (x >= 22 & x < 24) {
      return(0.6)
    } else if (x >= 24 & x < 26) {
      return(0.1)
    } else {
      return(0)
    }
  }
  
  custom_func_list <- replicate(length(trial_stim), custom_post_func, simplify = FALSE)
  custom_start <- runif(1, 20, 25)
  
  abs_model <- Zhu23ABS$new(width = 1, n_chains = 3, nd_time = 0.3, s_nd_time = 0.2, lambda = 10,
                            custom_distr = custom_func_list, custom_start = custom_start)
  abs_model$simulate(stopping_rule = 'fixed', n_sample = 5, trial_stim = trial_stim)
  expect_equal(abs_model$sim_results$samples[[1]][1], custom_start)
})


test_that("The custom starting point with relative stopping rule",{
  trial_stim <- factor(sample(c('left', 'right'), 2, TRUE))
  custom_post_left <- function(x){
    if (x >= -3 & x < -1){
      return(0.25 * x + 0.75)
    } else if (x >= -1 & x < 0) {
      return(-0.25 * x + 0.25)
    } else {
      return (0)
    }
  }
  
  custom_post_right <- function(x){
    if (x >= -1 & x < 1){
      return(0.25 * x + 0.25)
    } else if (x >= 1 & x < 3) {
      return(-0.25 * x + 0.75)
    } else {
      return (0)
    }
  }
  
  custom_func_list <- lapply(trial_stim, function(stim) ifelse(stim=='left', custom_post_left, custom_post_right))
  
  custom_start <- runif(1, -3, 3)
  
  abs_model2 <- Zhu23ABS$new(width=1, n_chains = 3, nd_time = 0.3, s_nd_time = 0.2, lambda = 10,
                             custom_distr = custom_func_list, custom_start = custom_start)
  abs_model2$simulate(stopping_rule = 'relative', delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim)
  expect_equal(abs_model2$sim_results$samples[[1]][1], custom_start)
})


test_that("starting points overwrites custom_start with fixed stopping rule", {
  trial_stim <- sample(20:25, 2, replace=TRUE)
  custom_post_func <- function(x){
    if (x >= 19 & x < 22){
      return(0.3)
    } else if (x >= 22 & x < 24) {
      return(0.6)
    } else if (x >= 24 & x < 26) {
      return(0.1)
    } else {
      return(0)
    }
  }
  
  custom_func_list <- replicate(length(trial_stim), custom_post_func, simplify = FALSE)
  custom_start <- runif(1, 20, 25)
  
  while (TRUE) {
    start_point <- runif(length(trial_stim), 20, 25)
    if (start_point[1] != custom_start){
      break
    }
  }
  
  abs_model <- Zhu23ABS$new(width = 1, n_chains = 3, nd_time = 0.3, s_nd_time = 0.2, lambda = 10,
                            custom_distr = custom_func_list, custom_start = custom_start)
  abs_model$simulate(stopping_rule = 'fixed', start_point = start_point, n_sample = 5, trial_stim = trial_stim)
  first_sample <- sapply(abs_model$sim_results$samples, function(samples) samples[1])
  expect_equal(first_sample, start_point)
})


test_that("starting points overwrites custom_start with fixed stopping rule",{
  trial_stim <- factor(sample(c('left', 'right'), 2, TRUE))
  custom_post_left <- function(x){
    if (x >= -3 & x < -1){
      return(0.25 * x + 0.75)
    } else if (x >= -1 & x < 0) {
      return(-0.25 * x + 0.25)
    } else {
      return (0)
    }
  }
  
  custom_post_right <- function(x){
    if (x >= -1 & x < 1){
      return(0.25 * x + 0.25)
    } else if (x >= 1 & x < 3) {
      return(-0.25 * x + 0.75)
    } else {
      return (0)
    }
  }
  
  custom_func_list <- lapply(trial_stim, function(stim) ifelse(stim=='left', custom_post_left, custom_post_right))
  custom_start <- runif(1, -3, 3)
  
  while (TRUE) {
    start_point <- runif(length(trial_stim), 20, 25)
    if (start_point[1] != custom_start){
      break
    }
  }
  
  abs_model2 <- Zhu23ABS$new(width=1, n_chains = 3, nd_time = 0.3, s_nd_time = 0.2, lambda = 10,
                             custom_distr = custom_func_list, custom_start = custom_start)
  abs_model2$simulate(stopping_rule = 'relative', start_point = start_point, delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim)
  first_sample <- sapply(abs_model2$sim_results$samples, function(samples) samples[1])
  expect_equal(first_sample, start_point)
})