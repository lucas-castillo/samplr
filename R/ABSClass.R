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
    #' @field delta The stop rule. It should be an integer.
    delta = NULL,
    
    #' @field lambda The rate parameter of the gamma distribution for decision time.
    lambda = NULL,
    
    #' @field n_chains The number of chains of the sampler. It should be an integer.
    n_chains = NULL,
    
    #' @field width The proposal width of the sampler.
    width = NULL,
    
    #' @field distr_name The type of the posterior hypothesis distribution.
    distr_name = NULL,
    
    
    #' @description
    #' Create a new 'CoreABS' object.
    #'
    #' @param delta The threshold of the relative difference. It should be an integer.
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param width The proposal width of the sampler.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #'
    #' @return A new 'CoreABS' object.
    #'
    initialize = function(delta, lambda, n_chains, width, distr_name='norm'){
      # Check variable types
      stopifnot("delta should be an integer."=(delta%%1 == 0))
      stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
      stopifnot("n_chains should be an integer."=(n_chains%%1 == 0))
      stopifnot("width should be a single numeric value."=(is.numeric(width) && length(width) == 1))
      # add checks for `distr_name`
      self$delta <- delta
      self$lambda <- lambda
      self$n_chains <- n_chains
      self$width <- width
      self$distr_name <- distr_name
    },
    
    two_alt_force_choice = function(){
      stop("Method 'confidence' not implemented.", call. = FALSE)
    },
    
    estimate = function(){
      invisible(TRUE)
    },
    
    prob_judgment = function(){
      invisible(TRUE)
    }
)
)


#' @title Auto-correlated Bayesian Sampler by Zhu (2023)
#' 
#' @description
#' This Auto-correlated Bayesian Sampler (ABS) model is developed by Zhu.
#' 
#' @examples
#' abs_model <- Zhu23ABS$new(nd_time=0.3, s_nd_time=0.5, prior_on_resp=c(1, 0), delta=3, lambda=10, n_chains=5, width=1, distr_name='norm')
#'
Zhu23ABS <- R6::R6Class(
  "Zhu23ABS",
  inherit = CoreABS,
  public = list(
    
    #' @field nd_time The non-decision time.
    nd_time = NULL,
    
    #' @field s_nd_time The range of the non-decision time.
    s_nd_time = NULL,
    
    #' @field prior_on_resp The beta prior on responses
    prior_on_resp = NULL,
    
    #' @description
    #' Create a new 'Zhu23ABS' object.
    #' 
    #' @param nd_time The non-decision time.
    #' @param s_nd_time The range of the non-decision time. Default is 0, implying a fixed non-decision time.
    #' @param prior_on_resp The beta prior on responses
    #' @param delta The stopping rule of the sampling process. It should be an integer.
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param width The proposal width of the sampler.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #'
    #' @return A new 'Zhu23ABS' object.
    #'
    initialize = function(nd_time, s_nd_time=0, prior_on_resp=c(0, 0), delta, lambda, n_chains, width, distr_name='norm') {
      super$initialize(delta, lambda, n_chains, width, distr_name)
      
      stopifnot("nd_time should be a single numeric value."=(is.numeric(nd_time) && length(nd_time) == 1))
      stopifnot("s_nd_time should be a single numeric value."=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
      stopifnot("prior_on_resp should be a numeric vector with two values."=(is.numeric(prior_on_resp) && length(prior_on_resp) == 2))
      self$nd_time <- nd_time
      self$s_nd_time <- s_nd_time
      self$prior_on_resp <- prior_on_resp
    },
    
    
    #' @description
    #' This function is for simulating two-alternative-force choice tasks by Zhu23ABS.
    #' 
    #' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
    #' @param discrim The stimuli discriminability.
    #' @param trial_stim The feedback of each trial. It should be a numeric vector only consist of 0 (below the decision boundary) and/or 1 (above the decision boundary).
    #' @param mc3_iterations The number of iterations for each chunk. See details for more information.
    #' 
    #' @return A data frame with ten columns:
    #' \enumerate{
    #'  \item{trial: the index of trials;}
    #'  \item{samples: the samples of ABS sampler for the trial;}
    #'  \item{support: the relative difference between the numbers of supporting samples for each response;}
    #'  \item{length: the length of the sampling sequence;}
    #'  \item{response: the response predicted by ABS. 0 represents the response coded as below the decision boundary and 1 represent the opposite;}
    #'  \item{stimulus: the stimuli of the experiment. 0 represents the response coded as below the decision boundary and 1 represent the opposite;}
    #'  \item{accuracy: whether the response is the same as the feedback. 0 represents error, and 1 represents correct;}
    #'  \item{nd_time: the non-decision time for this trial;}
    #'  \item{rt: the response time, including both the non-decision and the decision time;}
    #'  \item{confidence: the confidence of the response.}
    #' }
    #' 
    #' @importFrom magrittr %>%
    #' 
    #' @examples
    #' abs_model <- Zhu23ABS$new(nd_time=0.3, s_nd_time=0.5, prior_on_resp=c(1, 0), delta=3, lambda=10, n_chains=5, width=1, distr_name='norm')
    #' tafc_sim <- abs_model$two_alt_force_choice(dec_bdry=0, discrim=1, trial_stim=c(0, 0, 1, 0, 1))
    #' 
    two_alt_force_choice = function(dec_bdry, discrim, trial_stim, mc3_iterations=100){
      #Check inputs
      stopifnot("dec_bdry should be a single numeric value."=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
      stopifnot("discrim should be a single numeric value."=(is.numeric(discrim) && length(discrim) == 1))
      stopifnot("trial_stim should be a vector with either 0 and/or 1."=all(trial_stim %in% c(0, 1)))
      stopifnot("mc3_iterations should be a single numeric value."=(is.numeric(mc3_iterations) && length(mc3_iterations) == 1))
      
      start_point <- stats::runif(self$n_chains, min=-3, max=3) %>%
        as.matrix()
      
      tafc_sim <- Zhu23ABS_tafc_cpp(
        start_point = start_point,
        trial_stim = trial_stim, 
        distr_name = self$distr_name, 
        mc3_iterations = mc3_iterations, 
        n_chains = self$n_chains, 
        dec_bdry = dec_bdry, 
        discrim = discrim,
        prior_on_resp = self$prior_on_resp,
        delta = self$delta, 
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        er_lambda = self$lambda,
        proposal_width = self$width
      )
      tafc_df <- data.frame(do.call(rbind, tafc_sim))
      return(tafc_df)
    },
    
    estimate = function(){
      invisible(TRUE)
    },
    
    
    #' @examples
    #' abs_model <- Zhu23ABS$new(nd_time=0.3, s_nd_time=0.5, prior_on_resp=c(1, 0), delta=3, lambda=10, n_chains=5, width=1, distr_name='norm')
    #' pj_sim <- abs_model$prob_judgment(trial_bdry = c(3, 3, 4, 4), trial_stim = c(1, 2, 3, 4))
    #' 
    prob_judgment = function(trial_bdry, trial_stim){
      #Check inputs
      stopifnot("trial_bdry should be a numeric vector."=(is.numeric(trial_bdry)))
      stopifnot("trial_stim should be a numeric vector."=(is.numeric(trial_stim)))
      stopifnot("The length of trial_bdry and trial_stim should be equal."=(length(trial_bdry) == length(trial_stim)))
      
      start_point <- stats::runif(self$n_chains, min=-3, max=3) %>%
        as.matrix()
      
      pj_sim = Zhu23ABS_pj_cpp(
        start_point = start_point,
        trial_stim = trial_stim, 
        distr_name = self$distr_name,
        n_chains = self$n_chains, 
        trial_bdry = trial_bdry,
        prior_on_resp = self$prior_on_resp,
        delta = self$delta, 
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        er_lambda = self$lambda,
        proposal_width = self$width
      )
      pj_df <- data.frame(do.call(rbind, pj_sim))
      return(pj_df)
    }
  )
)



#' @title Auto-correlated Bayesian Sampler by Li
#' 
#' @description
#' This Auto-correlated Bayesian Sampler (ABS) model is developed by Yun-Xiao Li following the ABS of Zhu.
#' 
#' @examples
#' 
#'
LiABS <- R6::R6Class(
  "LiABS",
  inherit = CoreABS,
  public = list(
    
    #' @field nd_time The non-decision time.
    nd_time = NULL,
    
    #' @field s_nd_time The range of the non-decision time.
    s_nd_time = NULL
  )
)



