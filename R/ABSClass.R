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
   
   #' @field n_chains The number of chains of the sampler. It should be an integer.
   n_chains = NULL,
   #' @field distr_name The type of the posterior hypothesis distribution.
   distr_name = NULL,
   #' @field start_point The start point of each trial
   start_point = NULL,
   #' @field nd_time The non-decision time.
   nd_time = NULL,
   #' @field s_nd_time The range of the non-decision time.
   s_nd_time = NULL,
   
   
   #' @description
   #' Create a new 'CoreABS' object.
   #'
   #' @param n_chains The number of chains of the sampler. It should be an integer.
   #' @param distr_name The type of the posterior hypothesis distribution.
   #' @param start_point The start point of each trial.
   #' @param nd_time The non-decision time.
   #' @param s_nd_time The range of the non-decision time. Default is 0, implying a fixed non-decision time.
   #'
   #' @return A new 'CoreABS' object.
   #'
   initialize = function(n_chains, distr_name='norm', start_point=NA, nd_time, s_nd_time){
     # Check variable types
     
     stopifnot("n_chains should be an integer."=(n_chains%%1 == 0))
     # add checks for `distr_name`
     stopifnot("nd_time should be a single numeric value."=(is.numeric(nd_time) && length(nd_time) == 1))
     stopifnot("s_nd_time should be a single numeric value."=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
     
     self$n_chains <- n_chains
     self$distr_name <- distr_name
     self$start_point <- start_point
     self$nd_time <- nd_time
     self$s_nd_time <- s_nd_time
   }
)
)


#' @title Auto-correlated Bayesian Sampler by Zhu (2023)
#' 
#' @description
#' This Auto-correlated Bayesian Sampler model \insertCite{@ABS, @zhuAutocorrelatedBayesian2023}{samplr} is developed by Zhu.
#' 
#'
#'@references
#'    \insertRef{zhuAutocorrelatedBayesian2023}{samplr}
#'
#' @export
#'
Zhu23ABS <- R6::R6Class(
  "Zhu23ABS",
  inherit = CoreABS,
  public = list(
    
    #' @field lambda The rate parameter of the gamma distribution for decision time.
    lambda = NULL,
    #' @field width The standard deviation of the proposal distribution for MC3
    width = NULL,
    
    
    #' @description
    #' Create a new 'Zhu23ABS' object.
    #' 
    #' @param nd_time The non-decision time.
    #' @param s_nd_time The range of the non-decision time. Default is 0, implying a fixed non-decision time.
    #' @param lambda The rate parameter of the gamma distribution for decision time.
    #' @param n_chains The number of chains of the sampler. It should be an integer.
    #' @param width The proposal width of the sampler.
    #' @param distr_name The type of the posterior hypothesis distribution.
    #' @param start_point The start point of each trial.
    #' 
    #' @return A new 'Zhu23ABS' object.
    #'
    #' @examples
    #' zhuabs <- Zhu23ABS$new(nd_time = 0.3, s_nd_time = 0.5, lambda = 10, n_chains = 5, width = 1, distr_name = 'norm')
    #' 
    initialize = function(nd_time, s_nd_time=0, lambda, n_chains, width, distr_name='norm', start_point=NA) {
      super$initialize(n_chains, distr_name, start_point, nd_time, s_nd_time)
      
      stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
      stopifnot("width should be a single numeric value."=(is.numeric(width) && length(width) == 1))
      
      self$lambda <- lambda
      self$width <- width
    },
    
    
    #' @description
    #' This function is for simulating two-alternative-force choice tasks by Zhu23ABS.
    #' 
    #' @param delta The relative difference between the number of samples supporting each hypothesis.
    #' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
    #' @param prior_on_resp The beta prior on responses. Default setting is Beta(1,1)
    #' @param discrim The stimuli discriminability.
    #' @param trial_stim The stimulus of each trial. It should be a factor only consisting of two levels: below and above the decision boundary.
    #' @param stim_depend The boolean variable that control whether the prior on responses changes regarding the last stimulus.
    #' @param max_iterations The maximum length of the MCREC sampler. The program will stop the sampling process after the length of the sampling sequence reaches to this limitation.
    #' 
    #' @return A data frame with seven columns:
    #' \enumerate{
    #'  \item{trial: The index of trials;}
    #'  \item{samples: The samples of ABS sampler for the trial;}
    #'  \item{response: The response predicted by ABS;}
    #'  \item{stimulus: The stimuli of the experiment;}
    #'  \item{accuracy: Whether the response is the same as the feedback. 0 represents error, and 1 represents correct;}
    #'  \item{rt: The response time, including both the non-decision and the decision time;}
    #'  \item{confidence: The confidence of the response.}
    #' }
    #' 
    #' 
    #' 
    #' @examples
    #' trial_stim <- factor(sample(c('left', 'right'), 5, TRUE))
    #' tafc_sim <- zhuabs$two_alt_force_choice(delta = 4, dec_bdry = 0, discrim = 1, trial_stim = trial_stim)
    #' 
    two_alt_force_choice = function(delta, dec_bdry, prior_on_resp = c(1,1), discrim, trial_stim, stim_depend=TRUE, max_iterations=1000){
      #Check inputs
      stopifnot("delta should be an integer."=(delta %% 1==0))
      stopifnot("The value of delta should be larger than the absolute difference within prior_on_resp."=(delta>abs(prior_on_resp[0] - prior_on_resp[1])))
      stopifnot("prior_on_resp should be a numeric vector with two values."=(is.numeric(prior_on_resp) && length(prior_on_resp) == 2))
      stopifnot("dec_bdry should be a single numeric value."=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
      stopifnot("discrim should be a single numeric value."=(is.numeric(discrim) && length(discrim) == 1))
      stopifnot("trial_stim should be a factor."=is.factor(trial_stim))
      stopifnot("trial_depend should be a boolean variable."=(isTRUE(stim_depend) || isFALSE(stim_depend)))
      stopifnot("max_iterations should be a single numeric value."=(is.numeric(max_iterations) && length(max_iterations) == 1))
      
      trial_stim_num <- as.numeric(trial_stim)
      stim_levels <- levels(trial_stim)
      
      if (any(!is.na(self$start_point))){
        stopifnot("start_point should be a numeric vector" = (is.numeric(self$start_point)))
        stopifnot("The length of start_point should equal to the length of the stimuli." = (length(self$start_point) == length(trial_stim)))
      }
      
      tafc_sim <- Zhu23ABS_cpp(
        task_id = 2,
        trial_stim = trial_stim, 
        distr_name = self$distr_name,
        proposal_width = self$width,
        n_chains = self$n_chains,
        provided_start_point = self$start_point,
        prior_on_resp = prior_on_resp,
        stop_rule = delta,
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        lambda = self$lambda,
        stim_depend = stim_depend,
        mc3_iterations = max_iterations,
        dec_bdry = dec_bdry, 
        discrim = discrim
      )
      tafc_df <- data.frame(do.call(rbind, tafc_sim))
      tafc_df$stimulus <- trial_stim
      tafc_df$response <- stim_levels[as.numeric(tafc_df$response)]
      tafc_df$accuracy <- as.numeric(tafc_df$accuracy)
      tafc_df$rt <- as.numeric(tafc_df$rt)
      tafc_df$confidence <- as.numeric(tafc_df$confidence)
      
      return(tafc_df)
    },
    
    
    #' @description
    #' This function is for simulating the estimate and confidence interval tasks by Zhu23ABS.
    #' 
    #' @param n_sample The fixed number of samples for each trial.
    #' @param trial_stim The stimulus of each trial.
    #' @param conf_level The required confidence level.
    #' 
    #' @return A data frame with ten columns:
    #' \enumerate{
    #'  \item{trial: The index of trials;}
    #'  \item{samples: The samples of ABS sampler for the trial;}
    #'  \item{stimulus: The stimuli of the experiment;}
    #'  \item{nd_time: The non-decision time for this trial;}
    #'  \item{conf_interval_l: The lower bound of the confidence interval with the given level;}
    #'  \item{conf_interval_u: The upper bound of the confidence interval with the given level;}
    #' }
    #' 
    #' @examples
    #' trial_stim <- c(25, 26, 10, 30)
    #' est_sim <- zhuabs$estimate(n_sample = 5, trial_stim = trial_stim, conf_level = 0.5)
    #'
    estimate = function(n_sample, trial_stim, conf_level=FALSE){
      
      sss_df <- private$single_stim_simulation(n_sample, trial_stim)
      sss_df$point_est <- sapply(sss_df$samples, function(samples) samples[length(samples)])
      
      if (conf_level){
        stopifnot("conf_level should be a single value between 0 and 1."=(is.numeric(conf_level) & length(conf_level) == 1 & conf_level >= 0 & conf_level <= 1))
        conf_interval <- t(sapply(sss_df$samples, function(samples) quantile(samples, probs = c((1-conf_level)/2, (1+conf_level)/2))))
        sss_df$conf_interval_l <- conf_interval[,1]
        sss_df$conf_interval_u <- conf_interval[,2]
      }
      return(sss_df)
    }
  ),
  
  private = list(
    single_stim_simulation = function(n_sample, trial_stim){
      #Check inputs
      stopifnot("trial_stim should be a numeric vector."=(is.numeric(trial_stim)))
      if (any(!is.na(self$start_point))){
        stopifnot("start_point should be a numeric vector" = (is.numeric(self$start_point)))
        stopifnot("The length of start_point should equal to the length of the stimuli." = (length(self$start_point) == length(trial_stim)))
      }
      
      sss_sim = Zhu23ABS_cpp(
        task_id = 1,
        trial_stim = trial_stim, 
        distr_name = self$distr_name,
        n_chains = self$n_chains,
        proposal_width = self$width,
        provided_start_point = self$start_point,
        stop_rule = n_sample,
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        lambda = self$lambda,
        prior_on_resp = c(1,1) # a place holder
      )
      sss_df <- data.frame(do.call(rbind, sss_sim))
      return(sss_df)
    }
  )
)