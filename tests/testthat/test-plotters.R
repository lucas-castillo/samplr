sequence <- c(
  0, -1, -2, 1, 0, -6, -1, -2, -12, -10, -3, -5, -12, -18, 1, -3, -16,
  -5, -4, 7, -2, -17, -5, 7, -7, 9, -4, -24, 4, -6, 1, -7, -16, 0, 6, 
  -7, -23, -24, 17, 0, 10, -10, -9, 12, 14, -15, -31, -12, 11, -26, 
  -18, 16, 7, 10, -15, -15, -27, 8, 11, 10, 19, -28, -22
)

test_that("Levy Flights", {
  res <- plot_levy(sequence,F)
  
  # It's a list of length 4 with named outputs
  expect_type(res, "list")
  expect_equal(length(res), 4)
  expect_equal(names(res), c("fx", "fy", "slope", "coef"))
  
  # Check results
  expect_equal(res[["coef"]], c(-0.4249660, 0.6237538))
  
  ## Test that it works the same with matrices
  res2 <- plot_levy(matrix(sequence), F)
  expect_equal(res, res2)
})

test_that("PSD", {
  res <- plot_PSD(sequence,F)
  
  # It's a list of length 3 with named outputs
  expect_type(res, "list")
  expect_equal(length(res), 3)
  expect_equal(names(res), c("log_freq", "log_psd", "polyfit"))
  
  # Check results
  expect_equal(res[["polyfit"]], c(0.5096142, 2.4102135))
})

test_that("Sigma Scaling", {
  ## Error if bad input
  expect_error(plot_sigma_scaling(matrix(sequence, ncol=3), F), "Please input a one-dimensional vector")
  
  res <- plot_sigma_scaling(sequence, F)

  # It's a vector
  expect_type(res, "double")
  expect_equal(length(res), round(length(sequence) / 10))

  # Check results
  expect_equal(res, c(16.40681, 20.13771, 19.16075, 17.64051, 15.55798, 17.40862), tolerance = .00001)
})

