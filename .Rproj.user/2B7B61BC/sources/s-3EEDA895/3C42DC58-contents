checkDistrInfo <- function(distr_name, distr_params, start){
  no_name = is.null(distr_name);
  no_params = is.null(distr_params)

  if (no_name && no_params){
    stop("No distribution provided. Please provide a value for distr_name and for distr_params")
  } else if (no_name){
    stop("No distribution name provided. Please provide a value for distr_name.")
  } else if (no_params){
    stop("No distribution parameters provided. Please provide a value for distr_params")
  }

  # if info provided, make sense of it
  # possible distr names
  names_cont <- c("unif", "norm","lnorm", "gamma", "beta", "nbeta", "chisq", "nchisq", "t", "nt", "f", "nf", "cauchy", "exp", "logis", "weibull",
                  "4beta", "lst", "truncnorm", "trunct", "trunclst", "triangular")

  names_cont_mv <- c("mvnorm", "mvt")

  names_discr <- c("binom", "nbinom", "nbinom_mu", "pois", "geom", "hyper", "wilcox", "signrank")

  # possible distr parameters
  parameters_cont <- c(2, 2, 2, 2, 2, 3, 1, 2, 1, 2, 2, 3, 2, 1, 2, 2,
                       4, 3, 4, 3, 5, 3)
  parameters_discr <- c(2,2,2,1,1,3,2,1)
  parameters_cont_mv <- c(2, 3)

  isContinuous = is.element(distr_name, names_cont) || is.element(distr_name, names_cont_mv);
  isDiscrete = is.element(distr_name, names_discr);

  isValidName = isContinuous || isDiscrete

  isValidParameters = FALSE

  if (isValidName){
    if (isContinuous){
      if (is.element(distr_name, names_cont)){
        isValidParameters = parameters_cont[match(distr_name, names_cont)] == length(distr_params)
      } else{
        isValidParameters = parameters_cont_mv[match(distr_name, names_cont_mv)] == length(distr_params)
      }
    } else{
      isValidParameters = parameters_discr[match(distr_name, names_discr)] == length(distr_params)
    }
  }


  if (!isValidName){
    stop("Distribution name given is not currently supported")
  } else if (!isValidParameters){
    stop("Parameters given do not match Distribution name given (either too many, or too few, or not given in a list?)")
  }
  if (!is.element(distr_name, names_cont_mv)){
    if (length(start) != 1) {
      stop("Start position has more variables than the distribution")
    }
  } else {
    if (length(start) != 2){
      stop("Start position should have 2 variables")
    }
  }
  return (isDiscrete)
}

check_sigma_prop <- function(sigma_prop, n_dim){
  # If no variance for the proposal distribution was given...
  if (is.null(sigma_prop)){
    if (length(start) == 1){
      warning("The variance of the proposal distribution was not given and defaulted to 1")
    } else{
      warning("The variance of the proposal distribution was not given and defaulted to the identity matrix")
    }
    # give it the arbitrary value of the identity matrix
    return(diag(n_dim))
  } else{
    if (!is.matrix(sigma_prop) && length(sigma_prop) == 1 && n_dim == 1){
      return(matrix(sigma_prop))
    }
    return(sigma_prop)
  }
}


#' Markov Chain Monte Carlo Sampler
#'
#'
#' This sampler navigates the proposal distribution following a random walk. At each step, it generates a new proposal from a proposal distribution (in this case a Gaussian centered at the current position) and chooses to accept it or reject it following the Metropolis-Hastings rule: it accepts it if the density of the posterior distribution at the proposed point is higher than at the current point. If the current position is denser, it still may accept the proposal with probability `proposal_density / current_density`.
#'
#' As mentioned, the proposal distribution is a Normal distribution. Its mean is the current position, and its variance is equal to the `sigma_prop` parameter, which defaults to the identity matrix if not specified.
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param sigma_prop Covariance matrix of the proposal distribution. If sampling in 1D space, it can be instead a number.
#' @param iterations Number of iterations of the sampler.
#' @return A list containing
#' \enumerate{
#'  \item{the history of visited places (a n x d matrix, n = iterations; d = dimensions)}
#'  \item{acceptance ratio - the proportions of proposals that were accepted (numeric)}
#' }
#' @export
#'
#' @examples
#'
#' # Sample from a normal distribution
#' metropolis_hastings <- sampler_mcmc(distr_name = "norm", distr_params = c(0,1), start = 1, sigma_prop = diag(1))
#'
#' # Not giving a sigma_prop issues a warning, but the sampler runs anyway with a default value
#' metropolis_hastings <- sampler_mcmc(distr_name = "norm", distr_params = c(0,1), start = 1)
#'
#'
sampler_mcmc<- function(distr_name, distr_params, start, sigma_prop = NULL, iterations = 1024L){
  isDiscrete = checkDistrInfo(distr_name, distr_params, start)
  sigma_prop = check_sigma_prop(sigma_prop, length(start))


  return (sampler_mcmc_cpp(start, sigma_prop, iterations, distr_name, distr_params, discreteValues = isDiscrete))
}

#' Metropolis-coupled MCMC sampler (MC3)
#'
#' This sampler is a variant of MCMC in which multiple parallel chains are run at different temperatures. The chains stochastically swap positions which allows the coldest chain to visit regions far from its starting point (unlike in MCMC). Because of this, an MC3 sampler can explore far-off regions, whereas an MCMC sampler may become stuck in a particular point of high density.
#'
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param nChains Number of chains to run.
#' @param sigma_prop Covariance matrix of the proposal distribution. If sampling in 1D space, it can be instead a number.
#' @param delta_T numeric, >1. Temperature increment parameter. The bigger this number, the steeper the increase in temperature between the cold chain and the next chain
#' @param swap_all Boolean. If true, every iteration attempts floor(nChains / 2) swaps. If false, only one swap per iteration.
#' @param iterations Number of iterations of the sampler.
#'
#' @export
#'
#' @examples
#'
#' # Sample from a normal distribution
#' mc_3 <- sampler_mc3(distr_name = "norm", distr_params = c(0,1), start = 1, sigma_prop = diag(1))
sampler_mc3<- function(distr_name, distr_params, start, nChains = 6, sigma_prop = NULL, delta_T = 4, swap_all = TRUE, iterations = 1024L){
  isDiscrete = checkDistrInfo(distr_name, distr_params, start)
  sigma_prop = check_sigma_prop(sigma_prop, length(start))



  samplerResults <- sampler_mc3_cpp(start, nChains, sigma_prop, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues = isDiscrete)

  M <- array(0, dim = (c(iterations, length(start), nChains)))
  for (i in 1:nChains){
    start = 1 + (i-1) * iterations
    end = start + iterations - 1
    M[,,i] = samplerResults[[1]][start:end,]
  }

  swap_hist = samplerResults[[3]][1:samplerResults[[4]][2], ]
  # colnames(swap_hist) <- c("Iteration", "Chain 1", "Chain 2")



  return(list(
    "Samples" = M,
    "Beta Values" = samplerResults[[2]],
    "Swap History" = swap_hist,
    "Swap Acceptance Ratio" = samplerResults[[4]][1] / samplerResults[[4]][2],
    "Acceptance Ratios" = samplerResults[[5]]
  ))

}

#' Hamiltonian Monte-Carlo Sampler.
#'
#' Hamiltonian Monte-Carlo, also called Hybrid Monte Carlo, is a sampling algorithm that uses Hamiltonian Dynamics to approximate a posterior distribution. Unlike MCMC and MC3, HMC uses not only the current position, but also a sense of momentum, to draw future samples. An introduction to HMC can be read [here](http://arxiv.org/abs/1701.02434)
#'
#'
#' This implementations assumes that the momentum is drawn from a normal distribution with mean 0 and identity covariance matrix (p ~ N (0, I) )
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param iterations Number of iterations of the sampler.
#'
#' @export
#' @examples
#'
#' HMC <- sampler_hmc(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)
sampler_hmc <- function(distr_name, distr_params, start, epsilon=.5, L=10, iterations=1024) {
  isDiscrete = checkDistrInfo(distr_name, distr_params, start)


  samplerResults <- sampler_hmc_cpp(start, distr_name, distr_params, epsilon, L, iterations)

  return(
    list(
      "Samples" = samplerResults[[1]],
      "Momentums" = samplerResults[[2]],
      "Acceptance Ratio" = samplerResults[[3]]
       )
    )

}

#' No U-Turn Sampler.
#'
#' Adapted from Hoffman and Gelman (2014). The No U-Turn Sampler (NUTS) aims to eliminate the need to set a number of steps L that is present in Hamiltonian Monte Carlo, which may lead to undesirable behaviour in HMC if not set correctly.NUTS does so by recursively building a set of candidate points that span the target distribution, and stopping when it starts to double back (hence its name). More information can be found [here](https://arxiv.org/abs/1111.4246)
#'
#'
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param delta_max Measure of the required accuracy of the simulation. The authors recommend a large value (1000)
#' @param iterations Number of iterations of the sampler.
#'
#' @export
#'
#' @examples
#' NUTS <- sampler_nuts(distr_name = "norm", distr_params = c(0,1), start = 1)
sampler_nuts <- function(distr_name, distr_params, start, epsilon=.5, delta_max=1000, iterations=1024) {
  isDiscrete = checkDistrInfo(distr_name, distr_params, start)

  samplerResults <- sampler_nuts_cpp(start, distr_name, distr_params, epsilon, iterations, delta_max)

  return(
    list(
      "Samples" = samplerResults[[1]]
    )
  )

}
