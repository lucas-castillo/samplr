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
    #' @field delta The threshold of the relative difference. It should be an integer.
    delta = NULL,
    
    #' @field lambda The rate parameter of the gamma distribution for decision time.
    lambda = NULL,
    
    #' @field n_chains The number of chains of the sampler. It should be an integer.
    n_chains = NULL,
    
    #' @field width The proposal width of the sampler.
    width = NULL,
    
    #' @field distr_name The type of the posterior hypothesis distribution.
    distr_name='norm',
    
    #' @field mc3_iterations The number of iterations for each chunk. See details for more information.
    mc3_iterations=100,
    
    
    #' @description
    #' Create a new 'CoreABS' object.
    #'
    #' @param delta The threshold of the relative difference. It should be an integer.
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param width The proposal width of the sampler.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #' @param mc3_iterations The number of iterations for each chunk. See details for more information.
    #'
    #' @return A new 'CoreABS' object.
    #'
    initialize = function(delta, lambda, n_chains, width, distr_name='norm', mc3_iterations=100){
      # Check variable types
      stopifnot("delta should be an integer."=(delta%%1 == 0))
      stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
      stopifnot("n_chains should be an integer."=(n_chains%%1 == 0))
      stopifnot("width should be a single numeric value."=(is.numeric(width) && length(width) == 1))
      stopifnot("mc3_iterations should be a single numeric value."=(is.numeric(mc3_iterations) && length(mc3_iterations) == 1))
      # add checks for `distr_name`
      self$delta <- delta
      self$lambda <- lambda
      self$n_chains <- n_chains
      self$width <- width
      self$distr_name <- distr_name
      self$mc3_iterations <- mc3_iterations
    },
    
    confidence = function() {
      stop("Method 'confidence' not implemented.", call. = FALSE)
    },
    
    two_alt_force_choice = function() {
      stop("Method 'twofc' not implemented.", call. = FALSE)
    }
)
)



#' @title Auto-correlated Bayesian Sampler by Li
#' 
#' @description
#' This Auto-correlated Bayesian Sampler (ABS) model is developed by Yun-Xiao Li following the ABS of Zhu. LiABS has an extra component of non-decision time when calculating the response time.
#' 
#' @examples
#' abs_model <- LiABS$new(delta=5, lambda=20, n_chains=3, width=1, nd_time = 0.1, s_nd_time = 0.3)
#' 
#' @export
#'
LiABS <- R6::R6Class(
  "LiABS",
  inherit = CoreABS,
  public = list(
    
    #' @field nd_time The non-decision time.
    nd_time = NULL,
    
    #' @field s_nd_time The range of the non-decision time.
    s_nd_time = NULL,
    
    #' @description
    #' Create a new 'LiABS' object.
    #' 
    #' @param nd_time The non-decision time.
    #' @param s_nd_time The range of the non-decision time.
    #' @param delta The threshold of the relative difference. It should be an integer.
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param width The proposal width of the sampler.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #' @param mc3_iterations The number of iterations for each chunk. See details for more information.
    #'
    #' @return A new 'LiABS' object.
    #'
    initialize = function(nd_time, s_nd_time, delta, lambda, n_chains, width, distr_name='norm', mc3_iterations=100) {
      super$initialize(delta, lambda, n_chains, width, distr_name, mc3_iterations)
      
      stopifnot("nd_time should be a single numeric value."=(is.numeric(nd_time) && length(nd_time) == 1))
      stopifnot("s_nd_time should be a single numeric value."=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
      self$nd_time <- nd_time
      self$s_nd_time <- s_nd_time
    },
    
    #' @description
    #' This function is for simulating two-alternative-force choice tasks by LiABS.
    #' 
    #' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
    #' @param discrim The stimuli discriminability.
    #' @param trial_fdbk The feedback of each trial. It should be a numeric vector only consist of 0 (below the decision boundary) and/or 1 (above the decision boundary).
    #' 
    #' @return A data frame with ten columns:
    #' \enumerate{
    #'  \item{trial: the index of trials;}
    #'  \item{samples: the samples of ABS sampler for the trial;}
    #'  \item{support: the relative difference between the numbers of supporting samples for each response;}
    #'  \item{length: the length of the sequence;}
    #'  \item{response: the response predicted by ABS. 0 represents the response coded as below the decision boundary and 1 represent the opposite;}
    #'  \item{feedback: the feedback or the true stimuli of the experiment. 0 represents the response coded as below the decision boundary and 1 represent the opposite;}
    #'  \item{accuracy: whether the response is the same as the feedback. 0 represents error, and 1 represents correct;}
    #'  \item{nd_time: the non-decision time for this trial;}
    #'  \item{rt: the response time, including both the non-decision and the decision time;}
    #'  \item{rept: whether the next trial repeats the same response as the current trial;}
    #' }
    #' 
    #' @importFrom magrittr %>%
    #' 
    #' @examples
    #' abs_model <- LiABS$new(delta=5, lambda=20, n_chains=3, width=1, nd_time = 0.1, s_nd_time = 0.3)
    #' tafc_sim <- abs_model$two_alt_force_choice(dec_bdry=0, discrim=1, trial_fdbk=c(0, 0, 1, 0, 1))
    #' 
    two_alt_force_choice = function(dec_bdry, discrim, trial_fdbk){
      #Check inputs
      stopifnot("dec_bdry should be a single numeric value."=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
      stopifnot("discrim should be a single numeric value."=(is.numeric(discrim) && length(discrim) == 1))
      stopifnot("trial_fdbk should be a vector with either 0 and/or 1."=all(trial_fdbk %in% c(0, 1)))
      
      start_point <- stats::runif(self$n_chains, min=-3, max=3) %>%
        as.matrix()
      
      abs_sim <- ABS_sampler_tafc_cpp(
        start_point = start_point,
        trial_fdbk = trial_fdbk, 
        distr_name = self$distr_name, 
        mc3_iterations = self$mc3_iterations, 
        n_chains = self$n_chains, 
        dec_bdry = dec_bdry, 
        d_sepn = discrim, 
        delta = self$delta, 
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        er_lambda = self$lambda,
        proposal_width = self$width
      )
      
      x <- abs_sim %>% 
        do.call(cbind, .) %>% 
        t %>% 
        data.frame() %>% 
        tibble::tibble() %>% 
        tidyr::unnest(-c(samples, support))
      
      if (length(trial_fdbk)>1){
        x$rept <- c(NA, x$response[2:(length(x$response))] == x$response[1:(length(x$response) - 1)])
      } else {
        x$rept <- NA
      }
      invisible(x)
    }
  )
)