# Z Identities ------------------------------------------------------------
test_that("Z Identities", {
  a=.1
  b=.1
  a_and_b=.1
  a_or_b=.1
  a_given_b=.1
  b_given_a=.1
  a_given_not_b=.1
  b_given_not_a=.1
  a_and_not_b=.1
  b_and_not_a=.1
  not_a=.1
  not_b=.1
  # probabilities below 0 
  expect_error(
    Z_identities(
      a=-20, b, a_and_b, a_or_b, a_given_b, b_given_a, 
      a_given_not_b, b_given_not_a, a_and_not_b, b_and_not_a, not_a, not_b
    )
  )
  # probabilities above 1 
  expect_error(
    Z_identities(
      a=20, b, a_and_b, a_or_b, a_given_b, b_given_a, 
      a_given_not_b, b_given_not_a, a_and_not_b, b_and_not_a, not_a, not_b
    )
    
  )
  # NAs returned if probabilities not given
  expect_true(
    any(is.na(    
      Z_identities(
        a=NULL, b, a_and_b, a_or_b, a_given_b, b_given_a, 
        a_given_not_b, b_given_not_a, a_and_not_b, b_and_not_a, not_a=NULL, not_b=NULL
      ))
    )
  )
})

