sequence <- c(
  0, -1, -2, 1, 0, -6, -1, -2, -12, -10, -3, -5, -12, -18, 1, -3, -16,
  -5, -4, 7, -2, -17, -5, 7, -7, 9, -4, -24, 4, -6, 1, -7, -16, 0, 6,
  -7, -23, -24, 17, 0, 10, -10, -9, 12, 14, -15, -31, -12, 11, -26,
  -18, 16, 7, 10, -15, -15, -27, 8, 11, 10, 19, -28, -22
)

test_that("Euclidean Distance", {
  expect_error(euc_d(matrix(1:3), 1:3), "The two items must be both vectors or matrices")
  expect_error(euc_d(1:4, 1:3), "The two points must have the same dimensions")
  
  expect_error(
    euc_d(matrix(1:10, ncol = 2), matrix(1:10, ncol = 5)), 
    "The two matrices must have the same number of columns, representing each of the dimensions of the points they contain"
  )

  expect_warning(
    euc_d(matrix(1:10, ncol = 2), matrix(1:20, ncol = 2)), 
    "Because the matrices have uneven number of rows, only the distances for the first 5 rows were done"
  )
  
  expect_equal(
    euc_d(matrix(rep(4, 10), ncol = 2), matrix(rep(5, 10), ncol = 2)), 
    expected = rep(sqrt(2), 5)
  )
  
})

test_that("Change 1D", {
  expect_warning(change_1d(1), "X must be longer than 1")
  expect_equal(change_1d(1:3), rep(1,2))
})
  

test_that("Levy Flights", {
  res <- calc_levy(sequence, F)

  # It's a list of length 4 with named outputs
  expect_type(res, "list")
  expect_equal(length(res), 4)
  expect_equal(names(res), c("fx", "fy", "slope", "coef"))

  # Check results
  expect_equal(res[["coef"]], c(-0.4249660, 0.6237538))

  ## Test that it works the same with matrices
  res2 <- calc_levy(matrix(sequence), F)
  expect_equal(res, res2)
  
  # Test that it works with mv sequences
  res3 <- calc_levy(matrix(rep(sequence, 2), ncol=2), F)
  expect_equal(res3[["coef"]], c(-0.4249660, 0.68771756))
  
  vdiffr::expect_doppelganger("Levy Plot", \(){calc_levy(sequence, T)})
})

test_that("PSD", {
  ## Error if bad input
  expect_error(calc_PSD(matrix(sequence, ncol = 3), F), "Please input a one-dimensional vector")

  res <- calc_PSD(sequence, F)

  # It's a list of length 3 with named outputs
  expect_type(res, "list")
  expect_equal(length(res), 3)
  expect_equal(names(res), c("log_freq", "log_psd", "polyfit"))

  # Check results
  expect_equal(res[["polyfit"]], c("slope"=0.5096142, "intercept"=2.4102135))
  vdiffr::expect_doppelganger("PSD Plot", \(){calc_PSD(sequence, T)})
})

test_that("Sigma Scaling", {
  ## Error if bad input
  expect_error(calc_sigma_scaling(matrix(sequence, ncol = 3), F), "Please input a one-dimensional vector")

  res <- calc_sigma_scaling(sequence, F)

  # It's a vector
  expect_type(res$sds, "double")
  expect_equal(length(res$sds), round(length(sequence) / 10))

  # Check results
  expect_equal(res$sds, c(16.40681, 20.13771, 19.16075, 17.64051, 15.55798, 17.40862), tolerance = .00001)
  
  vdiffr::expect_doppelganger("sig scaling", \(){calc_sigma_scaling(sequence, T)})
})

test_that("QQ Plotter", {
  ## Error if bad input
  expect_error(calc_qqplot(matrix(sequence, ncol = 3)), "Please input a one-dimensional vector")
  vdiffr::expect_doppelganger("qq Plot", \(){calc_qqplot(sequence, plot=T)})
  vdiffr::expect_doppelganger("qq change Plot", \(){calc_qqplot(sequence, change=F, plot=T)})
})

test_that("Autocorr Plotter", {
  ## Error if bad input
  expect_error(calc_autocorr(matrix(sequence, ncol = 3)), "Please input a one-dimensional vector")
  vdiffr::expect_doppelganger("change autocorrelation Plot", \(){calc_autocorr(sequence, plot=T)})
  vdiffr::expect_doppelganger("sequence autocorrelation Plot", \(){calc_autocorr(sequence, change=F, plot=T)})
  
})

test_that("Series Plotter", {
  ## Error if bad input
  expect_error(plot_series(matrix(sequence, ncol = 3)), "Please input a one-dimensional vector")
  
  vdiffr::expect_doppelganger("series plot", \(){plot_series(sequence)})
  vdiffr::expect_doppelganger("change plot", \(){plot_series(sequence, change = T)})
})

test_that("calc all", {
  vdiffr::expect_doppelganger("all plot", \(){calc_all(sequence, plot = T)})
})
