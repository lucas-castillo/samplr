# Euclidean Distance Calculator
#
# Calculates the distance between two points. If instead given two matrices A and B of size n x d, for n number of items of d dimensions each, it will return an n - 1 length vector with the distances between each of the items (i.e. the distance between the item A\[1,\] and B\[1,\]; between A\[2,\] and B\[2,\], and so on).
#
#
#
# It assumes that all dimensions have the same units or that all dimensions are normalized (mean = 0 and std = 1). Adapted from [here](https://hlab.stanford.edu/brian/euclidean_distance_in.html)
#
# @param x,y Vector or Matrix. If vector of length d, a d-dimensional point. If matrix of size n x d, n points of d dimensions.
#
# @return Numeric or Vector. Distance between points. If vector, distances between the ith points in the matrices.
#


# @examples
# # distance between two points
# euc_d(c(0,0), c(3,4))
#
# # ask for multiple distances
# M <- matrix(0, nrow = 3, ncol = 2)
# M2 <- matrix(0, nrow = 3, ncol = 2)
# for (i in 1:3){
#   M[i,] <- stats::runif(2) * 5
#   M2[i,] <- stats::runif(2) * 5
# }
#euc_d(M, M2)
euc_d <- function(x,y){
  euclidean_distance <- function(a,b){
    counter = 0
    for (i in 1:length(a)){
      counter <- counter + (a[i] - b[i]) ** 2
    }
    return(sqrt(counter))
  }

  if (!(is.matrix(x) & is.matrix(y)) & !(is.vector(x) & is.vector(y))){
    stop("The two items must be both vectors or matrices")
  }
  if (is.vector(x)){
    if (length(x) != length(y)){
      stop("The two points must have the same dimensions")
    } else {
      return(euclidean_distance(x,y))
    }
  } else {
    if (ncol(x) != ncol(y)){
      stop("The two matrices must have the same number of columns, representing each of the dimensions of the points they contain")
    } else{
      distances <- vector()
      for (i in 1:min(nrow(x),nrow(y))){
        distances[i] <- euclidean_distance(x[i,], y[i,])
      }
      if (nrow(x) != nrow(y)){
        warning(paste("Because the matrices have uneven number of rows, only the distances for the first",min(nrow(x),nrow(y)),"rows were done"))
      }
      return(distances)
    }
  }

  #

}

# Change Calculator
#
# Given a one-dimensional chain, it returns a 1d chain of step sizes and direction (e.g. given a chain \[1,4,1\], it returns \[3, -3\])
#
# @param X one-dimensional chain.
#
# @return one-dimensional chain of step sizes
#
# @examples
# chain <- stats::runif(10, 1, 10)
# change_1d(chain)
#
#
change_1d <- function(X){
  if (length(X) > 1){
    return(X[2:length(X)] - X[1:(length(X)-1)])
  } else{
    warning("X must be longer than 1")
  }
}


#' Levy Flights Calculator
#'
#' This function analyses if the length of the jumps the sampler is making (\eqn{l}) belongs to a Levy probability density distribution, \eqn{P(l) \approx l^{-\mu}}.
#'
#' Values of \eqn{\mu \approx 2} have been used to describe foraging in animals, and produce the most effective foraging \insertCite{viswanathan1999OptimizingSuccessRandom}{samplr}. See \insertCite{zhu2018MentalSamplingMultimodal;textual}{samplr} for a comparison of Levy Flight and PSD measures for different samplers in multimodal representations.
#'
#' @param chain Matrix of n x d dimensions, n = iterations, d = dimensions.
#' @param plot Boolean. plot Boolean. Whether to also plot the distance-frequency relationship.
#' 
#' @references
#'  \insertAllCited{}
#'
#' @return
#' If plot is true, it returns a simple plot with the log absolute difference in estimates and their frequency, as well as an estimate for the \eqn{\mu} parameter. If false it returns a list with what's required to make the plot.
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' calc_levy(chain1[[1]], plot=TRUE)
calc_levy <- function(chain, plot=FALSE){
  distances <- vector()

  if (is.vector(chain)){
    for (i in 2:length(chain)){
      distances[i-1] <- euc_d(chain[i], chain[i-1])
    }
  } else if (is.matrix(chain)){
    for (i in 2:nrow(chain)){
      distances[i-1] <- euc_d(chain[i,], chain[i-1,])
    }
  }#good
  nbins <- 70
  h <- graphics::hist(distances, plot = FALSE, breaks = seq(0, max(distances), length.out = nbins + 1))
  x <- h$breaks
  y <- h$counts

  x <- x[-length(y)] # equal length (as 1 more break than bins)
  # take only the breaks where the count is more than 0
  x <- x[y > 0]
  y <- y[y > 0]
  # ignore the first break
  y <- y[x > 0]
  x <- x[x > 0]
  # this allows to do log10
  x <- log(x, base = 10)
  y <- log(y, base = 10)



  fnBins <- 8;
  binWidths <- (max(x)-min(x))/(fnBins-1)
  fx <- pracma::zeros(fnBins,1)
  fy <- pracma::zeros(fnBins,1)
  up_bound <- x[1]
  k <- 0

  for (j in 1:(fnBins-1)){
    low_bound <- up_bound
    up_bound <- low_bound+binWidths
    avgY <- mean(y[x>=low_bound & x<up_bound])
    if (is.finite(avgY) & avgY>0){
      k <- k+1
      fx[k] <- low_bound
      fy[k] <- avgY
    }
  }

  fx <- fx[1:k]
  fy <- fy[1:k]
  coef <- pracma::polyfit(fx,fy,1)
  slope <- pracma::polyval(coef,fx)
  df <- data.frame(Fx = fx, Fy = fy, Slope = slope)
  if (plot) {
    caption <- latex2exp::TeX(paste("\u0024\\hat{\\mu}\u0024 =", round(-coef[1], 3)))
    x_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Absolute Difference in Estimates)")
    y_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Frequency)")
    plot(df$Fx, df$Fy, 
         main="Levy Flights", xlab=x_lbl, ylab=y_lbl, pch=19, sub=caption)
    abline(lm(df$Slope~df$Fx), col="blue", lwd=2)
  }    
  return(list(fx = fx, fy = fy, slope = slope, coef = coef))
}

#' Power Spectral Density Calculator
#'
#' This function estimates the log power spectral density against the log frequency, and calculates a slope \eqn{\alpha}.
#'
#' A number of studies have reported that cognitive activities contain a long-range slowly decaying autocorrelation. In the frequency domain, this is expressed as \eqn{S(f)} ~  \eqn{1/f^{-\alpha}}, with \eqn{f} being frequency, \eqn{S(f)} being spectral power, and \eqn{\alpha} \eqn{\epsilon} \eqn{[0.5,1.5]} is considered \eqn{1/f} scaling. See See \insertCite{zhu2018MentalSamplingMultimodal;textual}{samplr} for a comparison of Levy Flight and PSD measures for different samplers in multimodal representations.
#'
#'
#' @param chain Matrix of n x d dimensions, n = iterations, d = dimensions sequence
#' @param plot Boolean. Whether to return a plot or the elements used to make it.
#' @references 
#'  \insertAllCited{}
#' @return
#' Returns a list with log frequencies, log PSDs, and slope and intercept estimates.

#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' calc_PSD(chain1[[1]], plot= TRUE)
calc_PSD <- function(chain, plot = FALSE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  pd <- stats::spectrum(chain, plot = FALSE)
  lf <- log(pd$freq[pd$freq != 0], base = 10)
  lpsd <- log(pd$spec[pd$spec != 0], base = 10)
  windows <- 9
  lf_range <- (max(lf) - min(lf))*2/(windows+1)
  y <- pracma::zeros(windows,1)
  x <- pracma::zeros(windows,1)
  x[1] <- min(lf);
  y[1] <- mean(lpsd[lf>=x[1] & lf<x[1]+lf_range])
  x[1] <- x[1] + lf_range/2
  for (i in 1:(windows-1)){
    y[i+1] <- mean(lpsd[lf >= x[i] & lf < (x[i]+lf_range)])
    x[i+1] = x[i]+lf_range/2;
  }
  Fit <- pracma::polyfit(x, y, 1)

  if (plot) {
    df <- data.frame(lf = lf, lpsd = lpsd)


    caption <- latex2exp::TeX(paste("\u0024\\hat{\\alpha} = ", round(-Fit[1], 3)))
    x_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Frequency)")
    y_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(PSD)")
    plot(df$lf, df$lpsd, type="l",
         main="Power Spectral Density", xlab=x_lbl, ylab=y_lbl,sub=caption)
    abline(Fit[2], Fit[1], col="blue", lwd=2)
  } 
  names(Fit) <- c("slope", "intercept")
  return(list(
    log_freq = lf,
    log_psd = lpsd,
    polyfit = Fit
  ))
}

#' QQ-Plot Calculator
#'
#' Estimates values for a QQ plot of Empirical values against Theoretical values from a normal distribution, for either the chain points or the distances between successive points. Optionally, returns a plot as well as the values. 
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param change Boolean. If false, it calculates a qqplot of the given chain. If true, it creates a chain of step sizes.
#' @param plot Boolean. Whether to plot the QQ plot or just return the values. 

#' @return
#' A list with the theoretical and empirical quantiles, and the intercept and slope of the line connecting the points
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' calc_qqplot(chain1[[1]], plot = TRUE)
calc_qqplot <- function(chain, change = TRUE, plot=FALSE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }

  if (change){
    y <- change_1d(chain)
    title = "QQ Plot - Change"
  } else{
    y <- chain
    title = "QQ Plot - Values"
  }
  v <- qqnorm(y, plot.it = FALSE)
  probs <- c(.25, .75)
  x <- qnorm(probs)
  y <- as.vector(quantile(y, probs, names = FALSE, na.rm = TRUE))
  slope <- diff(y)/diff(x)
  int <- y[[1L]] - slope * x[[1L]]
  if (plot){
    plot(v$x, v$y, 
         xlab = "Theoretical Quantiles", ylab="Empirical Quantiles", main=title)
    abline(int, slope, col="blue", lwd=2)
  }
  list(
    "theoretical_quantiles"=v$x,
    "empirical_quantiles"=v$y,
    "intercept"=int,
    "slope"=slope
    
  )
}

#' Sigma Scaling Calculator
#'
#' Calculates the sigma scaling of the chain, and optionally plots the result. 
#'
# TODO: Ref needed
#' Sigma scaling is defined as the slope of the regression connecting log time lags and the standard deviation of value changes across time lags. Markets show values of 0.5.
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param plot Boolean. Whether to additionally plot the result.
#'
#' @return
#' A list containing the vector of possible lags, the sd of the distances at each lag, their log10 counterparts, and the calculated intercept and slope. 
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' calc_sigma_scaling(chain1[[1]], plot = TRUE)
#'
calc_sigma_scaling <- function(chain, plot=FALSE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  s_devs <- vector()
  maxLag = round(length(chain) / 10)
  for (i in 1:maxLag){
    distances <- (chain[1:(length(chain)-i)] - matrix(chain[-1:-i]))
    s_devs[i] <- stats::sd(stats::na.omit(distances))
  }

  logS <- log(s_devs, base = 10)
  logN <- log(1:maxLag, base = 10)
  model <- stats::lm(logS ~ logN)
  intercept <- model$coefficients[["(Intercept)"]]
  slope <- model$coefficients[["logN"]]
  
  
  if (plot){
    x_lims <- c(0, max(logN))
    y_lims <- c(min(logS), log10(maxLag^max(0.5,slope))+intercept)
    
    plot(logN, logS,
         xlim=x_lims, ylim=y_lims,
         xlab = latex2exp::TeX("\u0024log_{10} (\\Delta t)\u0024"),
         ylab= latex2exp::TeX("\u0024log_{10} (\\sigma(\\Delta x))\u0024"),
         main="Sigma Scaling", 
         sub = paste("slope =", round(slope, 3))
    )
    abline(intercept, slope, col="blue", lwd=2)  
  }
  
  return(list(
    "lag" = 1:maxLag,
    "log_lag" = logN,
    "sds" = s_devs,
    "log_sds" = logS,
    "intercept" = intercept,
    "slope" = slope
  ))  
}

#' Autocorrelation Calculator
#'
#' Calculates the autocorrelation of a given sequence, or of the size of the steps (returns).
#'
#' Markets display no significant autocorrelations in the returns of a given asset.
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param change Boolean. If true, plot the autocorrelation of the change series. If false, plot the autocorrelation of the given chain.
#' @param alpha Measure of Type I error - defaults to .05
#' @param lag.max Length of the x axis. How far to examine the lags.
#' @param plot Boolean. Whether to additionally plot the result.
#'
#' @return
#' A vector with the standard deviations at each lag
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' calc_autocorr(chain1[[1]], plot=TRUE)
calc_autocorr <- function(chain, change = TRUE, alpha = .05, lag.max = 100, plot=FALSE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  if (change){
    a <- stats::acf(chain[2:length(chain)] - chain[1:(length(chain)-1)], lag.max = lag.max, plot=FALSE)
  } else{
    a <- stats::acf(chain, lag.max = lag.max, plot=FALSE)
  }
  atc <- as.vector(a$acf)[-1]
  if (plot){
    plot(1:(min(lag.max, length(atc))), atc, 
         ylim=c(-1, 1), 
         xlab="Lag", ylab = "ACF", 
         main = paste("Autocorrelation of", ifelse(change, "Changes", "Values") )
         )
    if (!change){
      c_int <- 1 - alpha
      vcrit <- pracma::erfinv(c_int) * sqrt(2)
      lconf = -vcrit/sqrt(length(chain));
      upconf = vcrit/sqrt(length(chain));
      lines(1:lag.max, rep(lconf, lag.max), lty="dashed", col="blue")
      lines(1:lag.max, rep(upconf, lag.max), lty="dashed", col="blue")
    }
  }
  
  return(atc)
}

#' Series Plotter
#'
#' Plots the value of a one-dimensional series against the iteration where it occurred. Useful to see the general pattern of the chain (white noise, random walk, volatility clustering)
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param change Boolean. Whether to plot the series of values or the series of changes between values. 
#' @return
#' A series plot
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_series(chain1[[1]])
plot_series <- function(chain, change=FALSE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  if (change){
    y <- c(NA, diff(chain))
  } else{
    y <- chain
  }
  plot(1:length(chain), y,
       xlab = "Iteration", ylab = "Value",
       main = paste("Series of", ifelse(change, "Changes", "Values"))
       )
}
#' Diagnostics Wrapper
#'
#' Calculates all diagnostic functions in the samplr package for a given chain. Optionally, plots them. 
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param plot Boolean. Whether to additionally plot the diagnostics. 
#' @param acf.alpha,acf.lag.max Additional parameters to \link[samplr]{calc_autocorr}.
#' @return
#' A list with all diagnostic calculations (a list of lists); and optionally a grid of plots.
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' diagnostics <- calc_all(chain1[[1]])
#' names(diagnostics)
calc_all <- function(chain, plot=TRUE, acf.alpha=.05, acf.lag.max=100){
  # restore user parameters on function exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (plot){
    par(mfrow=c(3,3), mai=c(2,0.22,1.3,0.12), mar=c(5,3,3,2)+.1)
    plot_series(chain, change=FALSE)
    plot_series(chain, change=TRUE)
  }
  
  
  value_qqplot <- calc_qqplot(chain, change = FALSE, plot = plot)
  change_qqplot <- calc_qqplot(chain, change = TRUE, plot = plot)
  
  value_autocorr <- calc_autocorr(chain, change = FALSE, plot = plot, 
                                  alpha = acf.alpha, lag.max = acf.lag.max)
  change_autocorr <- calc_autocorr(chain, change = TRUE, plot = plot, 
                                   alpha = acf.alpha, lag.max = acf.lag.max)
  
  levy <- calc_levy(chain, plot = plot)
  PSD <- calc_PSD(chain, plot = plot)
  sigma_scaling <- calc_sigma_scaling(chain, plot = plot)
  
  return(list(
    "value_qqplot" = value_qqplot,
    "change_qqplot" = change_qqplot,
    "value_autocorr" = value_autocorr,
    "change_autocorr" = change_autocorr,
    "levy" = levy,
    "PSD" = PSD,
    "sigma_scaling" = sigma_scaling
  ))

}
