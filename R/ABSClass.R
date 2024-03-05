#' @title CoreABS Object
#' @description
#' This is the parent [R6][R6::R6Class] class of the Auto-correlated Bayesian Sampler \insertCite{@ABS, @zhuAutocorrelatedBayesian2023}{samplr}. It is a sequential sampling model assuming people draw autocorrelated samples from memory or beliefs, i.e., posterior of hypotheses.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @references
#'    \insertRef{zhuAutocorrelatedBayesian2023}{samplr}
#'
#' 
CoreABS <- R6::R6Class("CoreABS",
  public = list(
    
    #' @field lambda The rate parameter of the gamma distribution for decision time.
    lambda = NULL,
    #' @field n_chains The number of chains of the sampler. It should be an integer.
    n_chains = NULL,
    #' @field distr_name The type of the posterior hypothesis distribution.
    distr_name = NULL,
    
    
    #' @description
    #' Create a new 'CoreABS' object.
    #'
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #'
    #' @return A new 'CoreABS' object.
    #'
    initialize = function(lambda, n_chains, distr_name='norm'){
      # Check variable types
      stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
      stopifnot("n_chains should be an integer."=(n_chains%%1 == 0))

      # add checks for `distr_name`
      self$lambda <- lambda
      self$n_chains <- n_chains
      self$distr_name <- distr_name
    }
)
)
