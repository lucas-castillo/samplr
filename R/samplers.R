.checkNamesMatchParams <- function(distr_name, distr_params){

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
    stop(paste("Distribution", distr_name, "is not supported"))
  }



  if(isValidParameters){
    # if (length(start) == 1 && substr(returnString, 3, 4) == "mv"){
    #   stop("Start is length 1 but distribution is multivariate")
    # } else if (length(start) != 1 && substr(returnString, 3,4) == "uv"){
    #   stop("Distribution is univariate but start is not length 1")
    # }
    return(c(d_uv, c_mv)) ## logical vector c(isDiscrete, isMultivariate)
  } else{
    if (c_mv){
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_cont_mv[match(distr_name, names_cont_mv)], "parameters are expected in a list"))
    } else if (!d_uv){
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_cont[match(distr_name, names_cont)], "parameters are expected in a vector"))
    } else{
      stop(paste("Parameters given do not match distribution name given. For the", distr_name, "distribution,", parameters_discr[match(distr_name, names_discr)], "parameters are expected in a vector"))
    }
  }

}

.checkStart <- function(info, start_dim){
  if (start_dim == 1 && info[2]){
    stop("Start dimensions is 1 but distribution is multivariate")
  } else if (start_dim != 1 && !info[2]){
    stop("Distribution is univariate but start has more than one dimension")
  }
}

.checkWeights <- function(weights, distr_num){
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

.checkSigmaProp <- function(sigma_prop, n_dim){
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

.checkGivenInfo <- function(distr_name, distr_params, start, weights, caller, custom_density, sigma_prop = NULL){
  # returns list:
  # isDiscrete = distrInfo[[1]]
  # isMix = distrInfo[[2]]
  # weights = distrInfo[[3]]
  # sigma_prop = distrInfo[[4]]

  start_dim <- if (is.matrix(start)) ncol(start) else if (is.vector(start)) length(start) else NULL


  if (caller == "mh" || caller == "mc3"){
        sigma_prop = .checkSigmaProp(sigma_prop, start_dim)
  }
  if (is.null(custom_density)){
    if (is.null(distr_name) || is.null(distr_params)){
      stop("Please provide a distribution name and parameters")
    }
    ### Is mixture
    if (length(distr_name) > 1){
      # same amount of distr and parameters
      if (length(distr_name) != length(distr_params)){
        stop("The vector of distribution names and the list of distribution parameters must be of equal length")
      }

      for (i in 1:length(distr_name)){
        info = .checkNamesMatchParams(distr_name[i], distr_params[[i]])
        if (info[1]){
          stop(paste("Mixture Distributions are only supported if all distributions are continuous. Distribution", distr_name[[i]], "is discrete."))
        }
        if (caller != "grid"){
          .checkStart(info, start_dim)
        }
      }
      weights = .checkWeights(weights, length(distr_name))
      return(list(FALSE, TRUE, weights, sigma_prop))


    ### Is not mixture
    } else {
        info = .checkNamesMatchParams(distr_name, distr_params)
        if (caller != "grid"){
          .checkStart(info, start_dim)
        }
        weights = .checkWeights(weights, 1)
        return (list(info[1], FALSE, weights, sigma_prop))
    }
  } else{
    if (!is.null(distr_name) || !is.null(distr_params)){
      warning("Both (1) a custom density function and (2) distribution name and parameters were provided. Only the custom density function will be used.")
    }
    if (!is.null(weights)){
      warning("Weights were provided for a custom density function but will not be used.")
    }
    return(list(FALSE, FALSE, 1, sigma_prop))
  }


}

#' Density Plotter
#'
#' Plots a 2D map of the density of a distribution. If plot = FALSE, returns a dataframe with the density for each cell in the grid
#'
#' @param start Vector c(x, y) with the coordinates of the bottom-left corner of the map.
#' @param size Distance covered by the map. In other words, the top-right corner of the map has coordinates c(x + size, y + size)
#' @param cellsPerRow Number of cells to plot in every row. The higher, the more resolution
#' @param names Name of the distribution from which to sample from.
#' @param params Distribution parameters.
#' @param weights Distribution weights (if it's a mix of distributions)
#' @param customDensity Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @param plot Whether to return a plot or a dataframe with the density in each coordinate
#'
#' @return Density Plot or dataframe
#' @export
#'
#' @examples
#' # plot supported distribution
#' plot_2d_density(
#' c(-5, -5), 10, cellsPerRow = 100, names = c("mvnorm", "mvnorm"),
#' params = list(list(c(-2,1), diag(2)), list(c(2,1), diag(2)))
#' )
#'
#' # plot custom distribution
#' customDensity_r <- function(x){
#'     if (x[1] > 0 && x[1] < 3 && x[2] < 0 && x[2] > -3){
#'         return (1)
#'     } else {
#'         return (0)
#'     }
#'}
#' plot_2d_density(start = c(0,-4), size = 5, customDensity = customDensity_r)
#'
#'
plot_2d_density <- function(start, size, cellsPerRow = 50, names = NULL, params = NULL, weights = NULL, customDensity = NULL, plot = TRUE){

  xRange <- seq(from = start[1], to = start[1] + size, length.out = cellsPerRow)
  xxRange <- rep(xRange, cellsPerRow)

  yRange <- seq(from = start[2], to = start[2] + size, length.out = cellsPerRow)

  for (i in 1:cellsPerRow){
    if (i == 1){
      yyRange <- rep(yRange[i], cellsPerRow)
    } else {
      yyRange <- c(yyRange, rep(yRange[i], cellsPerRow))
    }
  }

  info <- .checkGivenInfo(names, params, NULL, weights, "grid", customDensity)
  weights = info[[3]]
  if (is.null(customDensity)){
    useCustom = FALSE
    customDensity = function(){}
    density <- gridDensity_cpp(names, params, length(names) > 1, weights, xxRange, yyRange, cellsPerRow, densityFunc = customDensity, useCustomDensity = useCustom)
  } else{
    useCustom = TRUE
    density <- gridDensity_cpp("names", c(0,1), FALSE, weights, xxRange, yyRange, cellsPerRow, densityFunc = customDensity, useCustomDensity = useCustom)
  }


  df <- data.frame(x = xxRange, y = yyRange, density = density)

  if (plot){
    map <- ggplot2::ggplot(df) +
      ggplot2::geom_raster(mapping = ggplot2::aes(x = x, y = y, fill = density)) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme_void()

    return(map)
  }
  else{
    return(df)
  }
}

#' Metropolis-Hastings (MH) Sampler
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
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#'
#' @return A list containing
#' \enumerate{
#'  \item{the history of visited places (an n x d matrix, n = iterations; d = dimensions)}
#'  \item{the history of proposed places (an n x d matrix, n = iterations; d = dimensions). Nothing is proposed in the first iteration (the first iteration is the start value) and so the first row is NA}
#'  \item{acceptance history - A logical vector where each item represents whether the proposal was accepted for that iteration or not. The first value is always FALSE, as the first iteration is the start point.}
#' }
#' @export
#'
#' @examples
#'
#' # Sample from a normal distribution
#' metropolis_hastings <- sampler_mh(
#'                        distr_name = "norm", distr_params = c(0,1),
#'                        start = 1, sigma_prop = diag(1)
#'                        )
#'
#' # Not giving a sigma_prop issues a warning, but the sampler runs anyway with a default value
#' metropolis_hastings <- sampler_mh(distr_name = "norm", distr_params = c(0,1), start = 1)
#'
#'

sampler_mh<- function(start, distr_name = NULL, distr_params = NULL, sigma_prop = NULL, iterations = 1024L, weights = NULL, custom_density = NULL){
  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "mh", custom_density, sigma_prop)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  sigma_prop = distrInfo[[4]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }
  res <- sampler_mh_cpp(start, sigma_prop, iterations, distr_name, distr_params, discreteValues = isDiscrete, isMix = isMix, weights = weights, custom_func = custom_density, useCustom = use_custom)
  res[[2]][1, ] <- NA

  return (res)
}

#' Metropolis-coupled MCMC sampler (MC3)
#'
#' This sampler is a variant of MH in which multiple parallel chains are run at different temperatures. The chains stochastically swap positions which allows the coldest chain to visit regions far from its starting point (unlike in MH). Because of this, an MC3 sampler can explore far-off regions, whereas an MH sampler may become stuck in a particular point of high density.
#'
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Either a vector or a matrix. If it is a vector, it will be the starting point of all the chains (with length = number of dimensions). If it's a matrix, every row will be the starting point of one chain (and so it must have as many rows as nChains, and as many columns as number of dimensions in the space).
#' @param nChains Number of chains to run.
#' @param sigma_prop Covariance matrix of the proposal distribution. If sampling in 1D space, it can be instead a number.
#' @param delta_T numeric, >1. Temperature increment parameter. The bigger this number, the steeper the increase in temperature between the cold chain and the next chain
#' @param swap_all Boolean. If true, every iteration attempts floor(nChains / 2) swaps. If false, only one swap per iteration.
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @export
#'
#' @examples
#'
#' # Sample from a normal distribution
#' mc_3 <- sampler_mc3(distr_name = "norm", distr_params = c(0,1), start = 1, sigma_prop = diag(1))
sampler_mc3<- function(start, distr_name = NULL, distr_params = NULL, sigma_prop = NULL, nChains = 6, delta_T = 4, swap_all = TRUE, iterations = 1024L, weights = NULL, custom_density = NULL){


  # Check nChains is integer
  if (floor(nChains) != nChains){
    stop("nChains must be a whole number")
  }
  # Chech start is matrix. If vector, make matrix.
  if (is.matrix(start)){
    if (nrow(start) != nChains){
      stop("The start matrix must have as many rows as nChains (and as many columns as number of dimensions in sampling space)")
    }
  } else{
    if (is.vector(start)){
      start = matrix(data = rep(start, nChains), nrow = nChains, byrow = T)
    }
  }
  start_dim <- if (is.matrix(start)) ncol(start) else if (is.vector(start)) length(start) else NULL
  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "mc3", custom_density, sigma_prop)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  sigma_prop = distrInfo[[4]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }

  samplerResults <- sampler_mc3_cpp(start, nChains, sigma_prop, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues = isDiscrete, isMix = isMix, weights = weights, custom_func = custom_density, useCustom = use_custom)

  M <- array(0, dim = (c(iterations, start_dim, nChains)))
  P <- array(0, dim = (c(iterations, start_dim, nChains)))

  for (i in 1:nChains){
    begin = 1 + (i-1) * iterations
    end = begin + iterations - 1

    M[,,i] = samplerResults[[1]][begin:end,]

    P[,,i] = samplerResults[[2]][begin:end,]
  }
  P[1, , ] <- NA

  swap_hist = samplerResults[[4]][1:samplerResults[[5]][2], ]
  if (nChains > 1){
    colnames(swap_hist) <- c("Iteration", "Chain 1", "Chain 2")
  }



  return(list(
    "Samples" = M,
    "Proposals" = P,
    "Beta Values" = samplerResults[[3]],
    "Swap History" = swap_hist,
    "Swap Acceptance Ratio" = samplerResults[[5]][1] / samplerResults[[5]][2],
    "Acceptance Ratios" = samplerResults[[6]]
  ))

}

#' Hamiltonian Monte-Carlo Sampler (HMC)
#'
#' Hamiltonian Monte-Carlo, also called Hybrid Monte Carlo, is a sampling algorithm that uses Hamiltonian Dynamics to approximate a posterior distribution. Unlike MH and MC3, HMC uses not only the current position, but also a sense of momentum, to draw future samples. An introduction to HMC can be read [here](http://arxiv.org/abs/1701.02434)
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
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @export
#' @examples
#'
#' HMC <- sampler_hmc(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)
sampler_hmc <- function(start, distr_name = NULL, distr_params = NULL, epsilon=.5, L=10, iterations=1024, weights = NULL, custom_density = NULL) {

  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "hmc", custom_density)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }




  if (isDiscrete){
    stop("Hamiltonian Monte-Carlo is not supported with discrete distributions")
  }

  samplerResults <- sampler_hmc_cpp(start, distr_name, distr_params, epsilon, L, iterations, isMix = isMix, weights = weights, custom_func = custom_density, useCustom = use_custom)

  return(
    list(
      "Samples" = samplerResults[[1]],
      "Momentums" = samplerResults[[2]],
      "Acceptance Ratio" = samplerResults[[3]]
       )
    )

}


#' Metropolis-Coupled Hamiltonian Monte Carlo (MCHMC)
#'
#' Metropolis-Coupled version of HMC, i.e. running multiple chains at different temperatures which stochastically swap positions.
#'
#' Metropolis-Coupled HMC does not support discrete distributions.
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @param nChains Number of chains to run.
#' @param delta_T numeric, >1. Temperature increment parameter. The bigger this number, the steeper the increase in temperature between the cold chain and the next chain
#' @param swap_all Boolean. If true, every iteration attempts floor(nChains / 2) swaps. If false, only one swap per iteration.


#' @export
#' @examples
#'
#' MCHMC <- sampler_mchmc(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)

sampler_mchmc <- function(start, distr_name = NULL, distr_params = NULL, epsilon=.5, L=10, nChains = 6, delta_T = 4, swap_all = TRUE, iterations = 1024L, weights = NULL, custom_density = NULL){

  # Check nChains is integer
  if (floor(nChains) != nChains){
    stop("nChains must be a whole number")
  }
  # # Chech start is matrix. If vector, make matrix.
  # if (is.matrix(start)){
  #   if (nrow(start) != nChains){
  #     stop("The start matrix must have as many rows as nChains (and as many columns as number of dimensions in sampling space)")
  #   }
  # } else{
  #   if (is.vector(start)){
  #     start = matrix(data = rep(start, nChains), nrow = nChains, byrow = T)
  #   }
  # }
  #
  # start_dim <- if (is.matrix(start)) ncol(start) else if (is.vector(start)) length(start) else NULL
  start_dim <- length(start)
  #
  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "hmc", custom_density)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }

  if (isDiscrete){
    stop("Metropolis-Coupled Recycled HMC is not supported with discrete distributions")
  }

  res <- sampler_mc_rec_cpp(
    start = start,
    nChains = nChains, # for non-tempered versions of HMC and HOR, set this to 1
    delta_T = delta_T, #
    swap_all = swap_all,
    iterations = iterations,
    distr_name = distr_name,
    distr_params=distr_params,
    discreteValues=isDiscrete,
    isMix=isMix,
    weights = weights,
    custom_func = custom_density,
    useCustom=use_custom,
    epsilon=epsilon,
    L=L,
    alpha=0 # for HMC, set this to zero
  )

  M <- array(0, dim = (c(iterations, start_dim, nChains)))
  P <- array(0, dim = (c(iterations, start_dim, nChains)))

  for (i in 1:nChains){
    begin = 1 + (i-1) * iterations
    end = begin + iterations - 1

    M[,,i] = res$chain[begin:end,]

    P[,,i] = res$momentums[begin:end,]
  }
  P[1, , ] <- NA



  swap_hist = res$swaps[1:res$swap_accepts_attempts[2], ]
  if (nChains > 1){
    colnames(swap_hist) <- c("Iteration", "Chain 1", "Chain 2")
  }

  return (
    list(
      "Samples" = M,
      "Momentums" = P,
      "Acceptance Ratios" = res$acceptances / (iterations - 1),
      "Beta Values" = res$beta,
      "Swap History" = swap_hist,
      "Swap Acceptance Ratio" = res$swap_accepts_attempts[1] / res$swap_accepts_attempts[2]
    )
  )
}



#' Recycled-Momentum HMC Sampler (REC)
#'
#' Recycled-Momentum HMC is a sampling algorithm that uses Hamiltonian Dynamics to approximate a posterior distribution. Unlike in standard HMC, proposals are autocorrelated, as the momentum of the current trajectory is not independent of the last trajectory, but is instead updated by a parameter alpha (see Details).
#'
#' While in HMC the momentum in each iteration is an independent draw,, here the momentum of the last utterance \eqn{p^{n-1}} is also involved. In each iteration, the momentum \eqn{p} is obtained as follows \deqn{p \gets \alpha \times p^{n-1} + (1 - \alpha^2)^{\frac{1}{2}} \times v}; where \eqn{v \sim N(0, I)}.
#'
#'
#' Recycled-Momentum HMC does not support discrete distributions.
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param alpha Recycling factor, from -1 to 1 (see Details).
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @export
#' @examples
#'
#' REC <- sampler_rec(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)

sampler_rec <- function(start, distr_name = NULL, distr_params = NULL, epsilon=.5, L=10, alpha=.1, iterations = 1024L, weights = NULL, custom_density = NULL){

  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "hmc", custom_density)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }

  if (isDiscrete){
    stop("Recycled HMC is not supported with discrete distributions")
  }

  res <- sampler_mc_rec_cpp(
    start = start,
    nChains = 1, # for non-tempered versions of HMC and HOR, set this to 1
    delta_T = 1, #
    swap_all = F,
    iterations = iterations,
    distr_name = distr_name,
    distr_params=distr_params,
    discreteValues=isDiscrete,
    isMix=isMix,
    weights = weights,
    custom_func = custom_density,
    useCustom=use_custom,
    epsilon=epsilon,
    L=L,
    alpha # for HMC, set this to zero
  )

  return (
    list(
      "Samples" = res$chain,
      "Momentums" = res$momentums,
      "Acceptance Ratio" = res$acceptances / iterations
    )
  )
}


#' Metropolis-Coupled Recycled-Momentum HMC Sampler (MCREC)
#'
#' Metropolis-Coupled version of Recycled-Momentum HMC, i.e. running multiple chains at different temperatures which stochastically swap positions.
#'
#' Metropolis-Coupled Recycled-Momentum HMC does not support discrete distributions.
#'
#' @param distr_name Name of the distribution from which to sample from.
#' @param distr_params Distribution parameters.
#' @param start Vector. Starting position of the sampler.
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param alpha Recycling factor, from -1 to 1 (see Details).
#' @param iterations Number of iterations of the sampler.
#' @param weights If using a mixture distribution, the weights given to each constituent distribution. If none given, it defaults to equal weights for all distributions.
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @param nChains Number of chains to run.
#' @param delta_T numeric, >1. Temperature increment parameter. The bigger this number, the steeper the increase in temperature between the cold chain and the next chain
#' @param swap_all Boolean. If true, every iteration attempts floor(nChains / 2) swaps. If false, only one swap per iteration.


#' @export
#' @examples
#'
#' MCREC <- sampler_mcrec(distr_name = "norm", distr_params = c(0,1), start = 1, epsilon = .01, L = 100)

sampler_mcrec <- function(start, distr_name = NULL, distr_params = NULL, epsilon=.5, L=10, alpha=.1, nChains = 6, delta_T = 4, swap_all = TRUE, iterations = 1024L, weights = NULL, custom_density = NULL){

  # Check nChains is integer
  if (floor(nChains) != nChains){
    stop("nChains must be a whole number")
  }
  # # Chech start is matrix. If vector, make matrix.
  # if (is.matrix(start)){
  #   if (nrow(start) != nChains){
  #     stop("The start matrix must have as many rows as nChains (and as many columns as number of dimensions in sampling space)")
  #   }
  # } else{
  #   if (is.vector(start)){
  #     start = matrix(data = rep(start, nChains), nrow = nChains, byrow = T)
  #   }
  # }
  #
  # start_dim <- if (is.matrix(start)) ncol(start) else if (is.vector(start)) length(start) else NULL
  start_dim <- length(start)
  #
  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "hmc", custom_density)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }

  if (isDiscrete){
    stop("Metropolis-Coupled Recycled HMC is not supported with discrete distributions")
  }

  res <- sampler_mc_rec_cpp(
    start = start,
    nChains = nChains, # for non-tempered versions of HMC and HOR, set this to 1
    delta_T = delta_T, #
    swap_all = swap_all,
    iterations = iterations,
    distr_name = distr_name,
    distr_params=distr_params,
    discreteValues=isDiscrete,
    isMix=isMix,
    weights = weights,
    custom_func = custom_density,
    useCustom=use_custom,
    epsilon=epsilon,
    L=L,
    alpha # for HMC, set this to zero
  )

  M <- array(0, dim = (c(iterations, start_dim, nChains)))
  P <- array(0, dim = (c(iterations, start_dim, nChains)))

  for (i in 1:nChains){
    begin = 1 + (i-1) * iterations
    end = begin + iterations - 1

    M[,,i] = res$chain[begin:end,]

    P[,,i] = res$momentums[begin:end,]
  }
  P[1, , ] <- NA



  swap_hist = res$swaps[1:res$swap_accepts_attempts[2], ]
  if (nChains > 1){
    colnames(swap_hist) <- c("Iteration", "Chain 1", "Chain 2")
  }

  return (
    list(
      "Samples" = M,
      "Momentums" = P,
      "Acceptance Ratios" = res$acceptances / (iterations - 1),
      "Beta Values" = res$beta,
      "Swap History" = swap_hist,
      "Swap Acceptance Ratio" = res$swap_accepts_attempts[1] / res$swap_accepts_attempts[2]
    )
  )
}

#' No U-Turn Sampler (NUTS)
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
#' @param custom_density Instead of providing names, params and weights, the user may prefer to provide a custom density function.
#' @export
#'
#' @examples
#' NUTS <- sampler_nuts(distr_name = "norm", distr_params = c(0,1), start = 1)
#'
sampler_nuts <- function(start, distr_name = NULL, distr_params = NULL, epsilon=.5, delta_max=1000, iterations=1024, weights = NULL, custom_density = NULL) {
  distrInfo = .checkGivenInfo(distr_name, distr_params, start, weights, "nuts", custom_density)
  isDiscrete = distrInfo[[1]]
  isMix = distrInfo[[2]]
  weights = distrInfo[[3]]
  if (is.null(custom_density)){
    custom_density <- function(x){}
    use_custom = FALSE
  } else{
    distr_name = ""
    distr_params = c(0,1)
    use_custom = TRUE
  }

  if (isDiscrete){
    stop("NUTS is not supported with discrete distributions")
  }
  samplerResults <- sampler_nuts_cpp(start, distr_name, distr_params, epsilon, iterations, delta_max, isMix = isMix, weights = weights, custom_func = custom_density, useCustom = use_custom)

  return(
    list(
      "Samples" = samplerResults[[1]]
    )
  )
}
