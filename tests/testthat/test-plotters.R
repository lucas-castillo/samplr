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

