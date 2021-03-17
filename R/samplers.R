checkNamesMatchParams <- function(distr_name, distr_params){

  names_cont <- c("unif", "norm","lnorm", "gamma", "beta", "nbeta", "chisq", "nchisq", "t", "nt", "f", "nf", "cauchy", "exp", "logis", "weibull",
                  "4beta", "lst", "truncnorm", "trunct", "trunclst", "triangular")

  names_cont_mv <- c("mvnorm", "mvt")

  names_discr <- c("binom", "nbinom", "nbinom_mu", "pois", "geom", "hyper", "wilcox", "signrank")

  # possible distr parameters
  parameters_cont <- c(2, 2, 2, 2, 2, 3, 1, 2, 1, 2, 2, 3, 2, 1, 2, 2,
                       4, 3, 4, 3, 5, 3)
  parameters_discr <- c(2,2,2,1,1,3,2,1)
  parameters_cont_mv <- c(2, 3)

  c_uv = is.element(distr_name, names_cont)
  c_mv = is.element(distr_name, names_cont_mv);
  d_uv = is.element(distr_name, names_discr);


  isValidParameters = FALSE

  # if name in one of these three, check parameters correct. Else stop and throw error.

  if (c_uv){
    isValidParameters = parameters_cont[match(distr_name, names_cont)] == length(distr_params)
  } else if(c_mv){
    isValidParameters = parameters_cont_mv[match(distr_name, names_cont_mv)] == length(distr_params) && is.list(distr_params)
  } else if (d_uv){
    isValidParameters = parameters_discr[match(distr_name, names_discr)] == length(distr_params)
  } else{
    stop(paste("Distribution", name, "is not supported"))
  }



  if(isValidParameters){
    # if (length(start) == 1 && substr(returnString, 3, 4) == "mv"){
    #   stop("Start is length 1 but distribution is multivariate")
    # } else if (length(start) != 1 && substr(returnString, 3,4) == "uv"){
    #   stop("Distribution is univariate but start is not length 1")
    # }
    return(c(d_uv, c_mv)) ## logical vector c(isDiscrete, isMultivariate)
  } else{
    if (substr(returnVector, 3, 4) == "mv"){
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_cont_mv[match(distr_name, names_cont_mv)], "parameters are expected in a list"))
    } else if (substr(returnVector, 1, 1) == "c"){
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_cont[match(distr_name, names_cont)], "parameters are expected in a vector"))
    } else{
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_discr[match(distr_name, names_discr)], "parameters are expected in a vector"))
    }
  }

}

checkStart <- function(info, start){
  dim = length(start)
  if (dim == 1 && info[2]){
    stop("Start is length 1 but distribution is multivariate")
  } else if (dim != 1 && !info[2]){
    stop("Distribution is univariate but start is longer than length 1")
  }
}

checkWeights <- function(weights, distr_num){
  if (distr_num != 1){
    if (is.null(weights)){
      weights = rep(1 / distr_num, length.out = distr_num)
      message("Equal weights given to all distributions")
    } else{
      if (distr_num != length(weights)){
        stop("The vector of distribution names and the vector list of distribution parameters must be of equal length")
      } else if (sum(weights) != 1){
        stop("The sum of the weights must equal 1")
      }
    }
  } else{
    if (!is.null(weights)){
      warning("Weights were provided for a single distribution and thus will be ignored")
    } else{
      weights = 1
    }
  }
  return(weights)
}

checkSigmaProp <- function(sigma_prop, n_dim){
  # If no variance for the proposal distribution was given...
  if (is.null(sigma_prop)){
    if (n_dim == 1){
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

checkGivenInfo <- function(distr_name, distr_params, start, weights, caller, sigma_prop = NULL){
  dim = length(start)
  if (caller == "mcmc" || caller == "mc3"){
    sigma_prop = checkSigmaProp(sigma_prop, dim)
  }

  ### Is mixture
  if (length(distr_name) > 1){
    # same amount of distr and parameters
    if (length(distr_name) != length(distr_params)){
      stop("The vector of distribution names and the list of distribution parameters must be of equal length")
    }

    for (i in 1:length(distr_name)){
      info = checkNamesMatchParams(distr_name[i], distr_params[[i]])
      if (info[1]){
        stop(paste("Mixture Distributions are only supported if all distributions are continuous. Distribution", distr_name[[i]], "is discrete."))
      }
      checkStart(info, start)
    }
    weights = checkWeights(weights, length(distr_name))
    return(list(FALSE, TRUE, weights, sigma_prop))


  ### Is not mixture
  } else {
      info = checkNamesMatchParams(distr_name, distr_params)
      checkStart(info, start)
      weights = checkWeights(weights, 1)
      return (list(info[1], FALSE, weights, sigma_prop))
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
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @param iterations Number of iterations of the sampler.
#'
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

sampler_mcmc<- function(distr_name, distr_params, start, sigma_prop = NULL, iterations = 1024L, weights = NULL){

  distrInfo = checkGivenInfo(distr_name, distr_params, start, weights, "mcmc", sigma_prop)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  sigma_prop = distrInfo[[4]]

  return (sampler_mcmc_cpp(start, sigma_prop, iterations, distr_name, distr_params, discreteValues = isDiscrete, isMix = isMix, weights = weights))
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
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#'
#' @export
#'
#' @examples
#'
#' # Sample from a normal distribution
#' mc_3 <- sampler_mc3(distr_name = "norm", distr_params = c(0,1), start = 1, sigma_prop = diag(1))
sampler_mc3<- function(distr_name, distr_params, start, nChains = 6, sigma_prop = NULL, delta_T = 4, swap_all = TRUE, iterations = 1024L, weights = NULL){
  distrInfo = checkGivenInfo(distr_name, distr_params, start, weights, "mc3", sigma_prop)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  sigma_prop = distrInfo[[4]]

  samplerResults <- sampler_mc3_cpp(start, nChains, sigma_prop, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues = isDiscrete, isMix = isMix, weights = weights)

  M <- array(0, dim = (c(iterations, length(start), nChains)))
  for (i in 1:nChains){
    start = 1 + (i-1) * iterations
    end = start + iterations - 1
    M[,,i] = samplerResults[[1]][start:end,]
  }

  swap_hist = samplerResults[[3]][1:samplerResults[[4]][2], ]
  colnames(swap_hist) <- c("Iteration", "Chain 1", "Chain 2")



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
#' This implementations assumes that the momentum is drawn from a normal distribution with mean 0 and identity covariance matrix (p ~ N (0, I)). Hamiltonian Monte Carlo does not support discrete distributions.
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @export
#' @examples
#'
#' HMC <- sampler_hmc(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)
sampler_hmc <- function(distr_name, distr_params, start, epsilon=.5, L=10, iterations=1024, weights = NULL) {
  distrInfo = checkGivenInfo(distr_name, distr_params, start, weights, "hmc")
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (isDiscrete){
    stop("Hamiltonian Monte-Carlo is not supported with discrete distributions")
  }

  samplerResults <- sampler_hmc_cpp(start, distr_name, distr_params, epsilon, L, iterations, isMix = isMix, weights = weights)

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
#' Like Hamiltonian MonteCarlo, it does not support discrete distributions.
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param delta_max Measure of the required accuracy of the simulation. The authors recommend a large value (1000)
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @export
#'
#' @examples
#' NUTS <- sampler_nuts(distr_name = "norm", distr_params = c(0,1), start = 1)
#'
sampler_nuts <- function(distr_name, distr_params, start, epsilon=.5, delta_max=1000, iterations=1024, weights = NULL) {
  distrInfo = checkGivenInfo(distr_name, distr_params, start, weights, "nuts")
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (isDiscrete){
    stop("NUTS is not supported with discrete distributions")
  }
  samplerResults <- sampler_nuts_cpp(start, distr_name, distr_params, epsilon, iterations, delta_max, isMix = isMix, weights = weights)

  return(
    list(
      "Samples" = samplerResults[[1]]
    )
  )
}
