
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

test_that("get_true_probabilities ", {
  expect_equal(
      sum(
        unlist(
          get_true_probabilities(
            .5, .2, .1, .3
            )
          )[1:4]
        )
  , 1)
})

test_that("Bayesian Sampler", {
  probs <- get_true_probabilities(.4, .4, .2, .2)
  res1 <- Bayesian_Sampler(
    a_and_b=probs$a_and_b,
    b_and_not_a = probs$b_and_not_a, 
    a_and_not_b = probs$a_and_not_b, 
    not_a_and_not_b = probs$not_a_and_not_b,
    beta = 2, 
    N = 20
  )
  res2 <- Bayesian_Sampler(
    a_and_b=probs$a_and_b,
    b_and_not_a = probs$b_and_not_a, 
    a_and_not_b = probs$a_and_not_b, 
    not_a_and_not_b = probs$not_a_and_not_b,
    beta = 2, 
    N = 200
  )
  # more samples = closer to true -- 
  expect_true(  abs(sum(unlist(res1[1:4])) - 1) > abs(sum(unlist(res2[1:4])) - 1))
  
  # Different vector lengths
  expect_error(
    Bayesian_Sampler(
      a_and_b=rep(probs$a_and_b, 2),
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = 2, 
      N = 200
    )
  )
  # Weird beta length
  expect_error(
    Bayesian_Sampler(
      a_and_b=probs$a_and_b,
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = c(2,3), 
      N = 200
    )
  )
  # Weird N length
  expect_error(
    Bayesian_Sampler(
      a_and_b=probs$a_and_b,
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = 2, 
      N = sample(1:10, 2)
    )
  )
  # Weird N2 length
  expect_error(
    Bayesian_Sampler(
      a_and_b=probs$a_and_b,
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = 2, 
      N = 200, N2 = c(100, 50)
    )
  )
  
  # N2 > N
  expect_error(
    Bayesian_Sampler(
      a_and_b=probs$a_and_b,
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = 2, 
      N = 200, N2 = 300
    )
  )
  
  # Probabilities don't add up to 1
  expect_error(
    Bayesian_Sampler(
      a_and_b=probs$a_and_b - .1,
      b_and_not_a = probs$b_and_not_a, 
      a_and_not_b = probs$a_and_not_b, 
      not_a_and_not_b = probs$not_a_and_not_b,
      beta = 2, 
      N = 200
    )
  )
  
  # You can do multiple trials in one call
  res <- matrix(unlist(  Bayesian_Sampler(
    a_and_b=rep(probs$a_and_b, 2),
    b_and_not_a = rep(probs$b_and_not_a, 2), 
    a_and_not_b = rep(probs$a_and_not_b, 2), 
    not_a_and_not_b = rep(probs$not_a_and_not_b, 2),
    beta = 2, 
    N = c(200, 2000)
  )), ncol=2, byrow=T)
  expect_true(abs(sum(res[1:4, 1]) - 1) > abs(sum(res[1:4, 2]) - 1))
  
  # N2 affects conjunctions but not base probabilities
  res <- Bayesian_Sampler(
    a_and_b=rep(probs$a_and_b, 2),
    b_and_not_a = rep(probs$b_and_not_a, 2), 
    a_and_not_b = rep(probs$a_and_not_b, 2), 
    not_a_and_not_b = rep(probs$not_a_and_not_b, 2),
    beta = 2, 
    N = c(200, 200), N2 = c(200, 5)
  )
  expect_true(sd(res$a) == 0 && sd(res$a_and_b) != 0)
  
  
  
})