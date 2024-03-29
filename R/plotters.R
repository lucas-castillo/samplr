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
    warning("The two items must be both vectors or matrices")
    return()
  }
  if (is.vector(x)){
    if (length(x) != length(y)){
      warning("The two points must have the same dimensions")
    } else {
      return(euclidean_distance(x,y))
    }
  } else {
    if (ncol(x) != ncol(y)){
      warning("The two matrices must have the same number of columns, representing each of the dimensions of the points they contain")
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


#' Levy Flights Plotter
#'
#' This plotter analyses if the length of the jumps the sampler is making (\eqn{l}) belongs to a Levy probability density distribution, \eqn{P(l) \approx l^{-\mu}}.
#'
#' Values of \eqn{\mu \approx 2} have been used to describe foraging in animals, and produce the most effective foraging [(Viswanathan et al., 1999)](https://www.nature.com/articles/44831). See [Zhu et al. 2018](https://dl.acm.org/doi/abs/10.5555/3327345.3327477) for a comparison of Levy Flight and PSD measures for different samplers in multimodal representations.
#'
#' @param chain Matrix of n x d dimensions, n = iterations, d = dimensions.
#' @param plot Boolean. plot Boolean. Whether to return a plot or the elements used to make it.
#'
#' @return
#' If plot is true, it returns a simple plot with the log absolute difference in estimates and their frequency, as well as an estimate for the \eqn{\mu} parameter. If false it returns a list with what's required to make the plot.
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_levy(chain1[[1]])
plot_levy <- function(chain, plot=TRUE){
  distances <- vector()

  if (is.vector(chain)){
    for (i in 2:length(chain)){
      distances[i-1] <- euc_d(chain[i], chain[i-1])
    }
  } else if (is.matrix(chain)){
    for (i in 2:nrow(chain)){
      distances[i-1] <- euc_d(chain[i], chain[i-1])
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
    caption <- latex2exp::TeX(paste("\u0024\\hat{\\mu}\u0024 =", -coef[1]))
    x_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Absolute Difference in Estimates)")
    y_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Frequency)")
    pl <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(Fx,Fy))
    pl <- pl +
      ggplot2::geom_point(data=df) +
      ggplot2::geom_path(mapping = ggplot2::aes(Fx, Slope), data=df, colour = "red") +
      ggplot2::labs(title = "Levy Flights", x = x_lbl, y = y_lbl, caption = caption)
    return(pl)
  } else{
    return(list(fx = fx, fy = fy, slope = slope, coef = coef))
  }
}

#' Power Spectral Density Plotter
#'
#' This function plots the log power spectral density against the log frequency, and calculates a slope \eqn{\alpha}.
#'
#' A number of studies have reported that cognitive activities contain a long-range slowly decaying autocorrelation. In the frequency domain, this is expressed as \eqn{S(f)} ~  \eqn{1/f^{-\alpha}}, with \eqn{f} being frequency, \eqn{S(f)} being spectral power, and \eqn{\alpha} \eqn{\epsilon} \eqn{[0.5,1.5]} is considered \eqn{1/f} scaling. See [Zhu et al. 2018](https://dl.acm.org/doi/abs/10.5555/3327345.3327477) for a comparison of Levy Flight and PSD measures for different samplers in multimodal representations.
#'
#'
#' @param chain Matrix of n x d dimensions, n = iterations, d = dimensions sequence
#' @param plot Boolean. Whether to return a plot or the elements used to make it.
#'
#' @return
#' If plot is, it returns a simple plot with the log PSD against the log frequency, as well as an estimate for the \eqn{\alpha} parameter. If false it returns a list with what's required to make the plot.

#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_PSD(chain1[[1]])
plot_PSD <- function(chain, plot = TRUE){
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


    caption <- latex2exp::TeX(paste("\u0024\\hat{\\alpha} = ", -Fit[1]))
    x_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(Frequency)")
    y_lbl <- latex2exp::TeX("\u0024log_{10} \u0024(PSD)")
    pl <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(lf, lpsd))
    pl <- pl + ggplot2::geom_line(mapping = ggplot2::aes(lf, lpsd), data = df) +
      ggplot2::geom_segment(mapping = ggplot2::aes(
        x = lf[1],
        y = pracma::polyval(Fit, lf)[1],
        xend =lf[length(lf)],
        yend = pracma::polyval(Fit, lf)[length(pracma::polyval(Fit, lf))]), colour = "red", size = 1) +
      ggplot2::labs(title = "Sigma Scaling", caption =caption, x = x_lbl, y = y_lbl)
    return(pl)

  } else{
    return(list(
      log_freq = lf,
      log_psd = lpsd,
      polyfit = Fit
    ))
  }
}

#' QQ-Plotter
#'
#' Plots a QQ plot of Empirical values against Theoretical values from a normal distribution. Can plot the chain points or the distances between successive points
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param change Boolean. If false, it plots a qqplot of the given chain. If true, it creates a chain of step sizes (using \link[SampleR]{change_1d})
#'
#' @return
#' QQ plot of Theoretical vs Empirical values
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_qqplot(chain1[[1]])
plot_qqplot <- function(chain, change = TRUE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }


  if (change){
    df <- data.frame(X = change_1d(chain))
    title = "QQ Plot - Change"
  } else{
    df <- data.frame(X = chain)
    title = "QQ Plot"
  }
  pl <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(sample = X))
  return(pl + ggplot2::stat_qq() + ggplot2::stat_qq_line(colour = "red") + ggplot2::ggtitle(title))
}

#' Sigma Scaling Plotter
#'
#' Plots a scaling of the sd in the distribution of price changes across time lags and returns the value of the slope
#'
#'Markets show sigma scaling exponents around 0.5.
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param plot Boolean. Whether to return a plot or the elements used to make it.
#'
#' @return
#' If plot is true, a sigma scaling plot. If false, a vector with the standard deviations at each lag
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_sigma_scaling(chain1[[1]])
#'
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_sigma_scaling(chain1[[1]], plot = FALSE)
plot_sigma_scaling <- function(chain, plot=TRUE){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  s_devs <- vector()
  maxLag = round(length(chain) / 10)
  for (i in 1:maxLag){
    distances <- (chain[1:(length(chain)-i)] - matrix(chain[-1:-i,]))
    s_devs[i] <- stats::sd(stats::na.omit(distances))
  }


  if (plot){
    logS <- log(s_devs, base = 10)
    logN <- log(1:maxLag, base = 10)
    sig_r <- pracma::polyfit(logN, logS, 1)
    x_lims <- c(0, max(logN))
    y_lims <- c(min(logS), log10(maxLag^max(0.5,sig_r[1]))+sig_r[2])

    df <- data.frame(logN = logN, logS = logS)

    pl <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(logN, logS))
    model <- stats::lm(logS ~ logN, data = df)
    intercept <- model$coefficients[["(Intercept)"]]
    slope <- model$coefficients[["logN"]]
    caption <- paste("Slope = ", slope)
    x_lbl <- latex2exp::TeX("\u0024log_{10} (\\Delta t)\u0024")
    y_lbl <- latex2exp::TeX("\u0024log_{10} (\\sigma(\\Delta t))\u0024")


    return(pl + ggplot2::geom_point() +
             ggplot2::geom_abline(ggplot2::aes(intercept = intercept, slope = slope), colour = "red") +
             ggplot2::labs(title = "Sigma Scaling", x = x_lbl, y = y_lbl, caption = caption) +
             ggplot2::xlim(x_lims) + ggplot2::ylim(y_lims))
  } else {
    return(s_devs)
  }

}

#' Autocorrelation Plotter
#'
#' Plots the autocorrelation of a given sequence, or of the size of the steps (returns).
#'
#' Markets display no significant autocorrelations in the returns of a given asset.
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param changeACF Boolean. If true, plot the autocorrelation of the change series. If false, plot the autocorrelation of the given chain.
#' @param alpha Measure of Type I error - defaults to .05
#' @param lag.max Length of the x axis. How far to examine the lags.
#'
#' @return
#' An autocorrelation plot
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_autocorr(chain1[[1]])
plot_autocorr <- function(chain, changeACF = TRUE, alpha = .05, lag.max = 100){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  if (changeACF){
    a <- stats::acf(chain[2:length(chain)] - chain[1:(length(chain)-1)], lag.max = lag.max, plot=FALSE)
    df <- data.frame(Lag = a$lag[-1], Autocorrelation = a$acf[-1])
    upperLine <- NULL
    lowerLine <- NULL
    title <- ggplot2::ggtitle("Change ACF")
  } else{
    c_int <- 1 - alpha
    a <- stats::acf(chain, lag.max = lag.max, plot=FALSE)
    df <- data.frame(Lag = a$lag, Autocorrelation = a$acf)
    vcrit <- pracma::erfinv(c_int) * sqrt(2)
    lconf = -vcrit/sqrt(length(chain));
    upconf = vcrit/sqrt(length(chain));
    upperLine <- ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = upconf), colour = "red", linetype = "dashed")
    lowerLine <- ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = lconf), colour = "red", linetype = "dashed")
    title <- ggplot2::ggtitle("Autocorrelation")
  }

  gpl <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(Lag, Autocorrelation))

  return (gpl + ggplot2::geom_line(data=df) + upperLine + lowerLine + ggplot2::ylim(c(-1,1)) + title)
  #   # ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = lconf), colour = "red", linetype = "dashed") +
  #   # ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = upconf), colour = "red", linetype = "dashed") +
  #     ggplot2::ylim(c(-1,1)) + ggplot2::ggtitle("Change ACF"))
}

#' Series Plotter
#'
#' Plots the value of a one-dimensional series against the iteration where it occurred. Useful to see the general pattern of the chain (white noise, random walk, volatility clustering)
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#'
#' @return
#' A series plot
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_series(chain1[[1]])
plot_series <- function(chain){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  df = data.frame(t = 1:length(chain), X = chain)
  return(ggplot2::ggplot(df, ggplot2::aes(t, X)) + ggplot2::geom_path(size=.1) + ggplot2::labs(title = "Series", x = "Iteration", y = "Value"))
}

#' Change Plotter
#'
#' Plots a change series against iterations. Useful to see if there is clustering of volatility in returns, like occurs in financial markets
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations

#' @return A plot of the change series
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_change(chain1[[1]])
plot_change <- function(chain){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  df = data.frame(t = 1:(length(chain)-1), X = change_1d(chain))
  return(ggplot2::ggplot(df, ggplot2::aes(t, X)) + ggplot2::geom_path(size=.1) + ggplot2::labs(title = "Change Series", x = "Iteration", y = "Change from previous"))
}

#' Plotter Wrapper
#'
#' Plots all the plot_* plots into a grid for ease.
#'
#' @param chain Vector of n length, where n is the number of trials or sampler iterations
#' @param title Title of the uberplot
#'
#' @return
#' A grid plotting all the plot_* functions
#' @export
#'
#' @examples
#' set.seed(1)
#' chain1 <- sampler_mh(1, "norm", c(0,1), diag(1))
#' plot_all(chain1[[1]])
plot_all <- function(chain, title = NULL){
  if (is.matrix(chain) && ncol(chain)>1){
    stop("Please input a one-dimensional vector")
  }
  p1 <- plot_series(chain)
  p2 <- plot_autocorr(chain)
  p3 <- plot_PSD(chain)
  p4 <- plot_sigma_scaling(chain)
  p5 <- plot_qqplot(chain)
  p6 <- plot_levy(chain)
  p7 <- plot_change(chain)

  full <- gridExtra::grid.arrange(p1,p2,p3,p4,p5, p6, p7, nrow = 3, ncol = 3, top = title, heights = rep(3,3), widths = rep(3,3))
  return(full)

}
