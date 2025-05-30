#' @title CoreABS Object
#' @description
#' This is the parent [R6][R6::R6Class] class of the Auto-correlated Bayesian Sampler \insertCite{@ABS, @zhu2024AutocorrelatedBayesianSampler}{samplr}. It is a sequential sampling model assuming people draw autocorrelated samples from memory or beliefs, i.e., posterior of hypotheses.
#' 
#' @importFrom R6 R6Class
#' @importFrom Rdpack reprompt
#' @keywords internal
#' @references
#'    \insertAllCited{}
#'
#' 
CoreABS <- R6::R6Class("CoreABS",
  public = list(
   
   #' @field n_chains an integer of the number of chains for the sampler.
   n_chains = NULL,
   #' @field nd_time a numeric value of the non-decision time (in seconds).
   nd_time = NULL,
   #' @field s_nd_time a numeric value of the inter-trial-variability of the non-decision time (in seconds).
   s_nd_time = NULL,
   #' @field distr_name a character string indicating the type of the posterior hypothesis distribution.
   distr_name = NULL,
   #' @field distr_params a numeric vector of the additional parameters for the posterior hypothesis distribution.
   distr_params = NULL,
   #' @field custom_distr a list of functions that define the posterior hypothesis distribution.
   custom_distr = NULL,
   #' @field custom_start a numeric value of the starting point if "custom_distr" is provided.
   custom_start = NULL,
   #' @field sim_results a data frame for saving the simulation results.
   sim_results = NULL,
  
   
   #' @description
   #' Create a new 'CoreABS' object.
   #' 
   #' @param n_chains an integer of the number of chains for the sampler.
   #' @param nd_time a numeric value of the non-decision time (in seconds).
   #' @param s_nd_time a numeric value of the inter-trial-variability of the non-decision time (in seconds).
   #' @param distr_name a character string indicating the type of the posterior hypothesis distribution. The package currently only supports `norm`, which represents normal distribution.
   #' @param distr_params a numeric vector of the additional parameters for the posterior hypothesis distribution.
   #' @param custom_distr a list of functions that define the posterior hypothesis distribution.
   #' @param custom_start a numeric value of the starting point if "custom_distr" is provided.
   #' 
   #' @return A new 'CoreABS' object.
   #'
   initialize = function(n_chains, nd_time, s_nd_time, distr_name = NULL, distr_params = NULL, custom_distr = NULL, custom_start = NULL){
     
     # Check variable types
     stopifnot('Argument "n_chains" should be a single integer.'=is.numeric(n_chains))
     stopifnot('Argument "n_chains" should be a single integer.'=((n_chains%%1 == 0) & length(n_chains) == 1))
     stopifnot('Argument "nd_time" should be a single numeric value.'=(is.numeric(nd_time) && length(nd_time) == 1))
     stopifnot('Argument "s_nd_time" should be a single numeric value.'=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
     
     # Check distributions
     if (is.null(distr_name) & is.null(custom_distr)) {
       stop('Both "distr_name" and "custom_distr" are empty.')
     }
     if (!is.null(distr_name)){
       stopifnot('Argument "distr_params" should be a numeric vector.'=(is.numeric(distr_params)))
     }
     if (!is.null(distr_name) & !is.null(custom_distr)){
       message('Both "distr_name" and "custom_distr" are provided. The distribution defined by "custom_distr" will be used as the posterior hypothesis.')
     }
     
     # Check `custom_start`
     if (!is.null(custom_distr)){
       stopifnot('Argument "custom_start" should be a numeric value.' = (is.numeric(custom_start) && length(custom_start) == 1))
     }
     
     self$n_chains <- n_chains
     self$nd_time <- nd_time
     self$s_nd_time <- s_nd_time
     self$distr_name <- distr_name
     self$distr_params <- distr_params
     self$custom_distr <- custom_distr
     self$custom_start <- custom_start
   }
),
  
  private = list(
    stopping_rule = NULL
  )
)


#' @title Auto-correlated Bayesian Sampler by Zhu (2023)
#' 
#' @description
#' This Auto-correlated Bayesian Sampler model \insertCite{@ABS, @zhu2024AutocorrelatedBayesianSampler}{samplr} is developed by Zhu.
#' 
#'
#' @references
#'    \insertAllCited{}
#'
#' @export
#'
Zhu23ABS <- R6::R6Class(
  "Zhu23ABS",
  inherit = CoreABS,
  public = list(
    
    #' @field width the standard deviation of the proposal distribution for MC3.
    width = NULL,
    #' @field lambda the rate parameter of the Erlang distribution for decision time.
    lambda = NULL,
    
    
    #' @description
    #' Create a new 'Zhu23ABS' object.
    #' 
    #' @param width a numeric value of the standard deviation of the proposal distribution for MC3.
    #' @param n_chains an integer of the number of chains for the sampler.
    #' @param nd_time a numeric value of the non-decision time (in seconds). When `s_nd_time` is not 0, `nd_time` represents the lower bound of the non-decision time.
    #' @param s_nd_time a numeric value of the inter-trial-variability of the non-decision time (in seconds).
    #' @param lambda a numeric value of the rate parameter of the Erlang distribution for decision time.
    #' @param distr_name a character string indicating the type of the posterior hypothesis distribution.
    #' @param distr_params a numeric vector of the additional parameters for the posterior hypothesis distribution.
    #' @param custom_distr a list of functions that define the posterior hypothesis distribution.
    #' @param custom_start a numeric value of the starting point if "custom_distr" is provided.
    #' 
    #' @return A new 'Zhu23ABS' object.
    #'
    #' @examples
    #' zhuabs <- Zhu23ABS$new(
    #'     width = 1, n_chains = 5, nd_time = 0.3, s_nd_time = 0.5, 
    #'     lambda = 10, distr_name = 'norm', distr_params = 1
    #' )
    #' 
    initialize = function(width, n_chains, nd_time, s_nd_time, lambda, distr_name = NULL, distr_params = NULL, custom_distr = NULL, custom_start = NULL) {
      super$initialize(n_chains, nd_time, s_nd_time, distr_name, distr_params, custom_distr, custom_start)
      
      stopifnot('Argument "lambda" should be a single numeric value.'=(is.numeric(lambda) && length(lambda) == 1))
      stopifnot('Argument "width" should be a single numeric value.'=(is.numeric(width) && length(width) == 1))
      
      self$lambda <- lambda
      self$width <- width
    },
    
    #' @description
    #' Simulate the ABS model.
    #' 
    #' 
    #' @param stopping_rule a character string indicating the stopping rule of ABS to be applied. Possible values are `"fixed"` and `"relative"`. See also `Details`.
    #' @param start_point a numeric vector setting the start point of each trial for the sampler. By default, it's set to `NA`, indicating that the starting point of the first trial is a random point from the posterior of hypotheses, and the starting points of subsequent trials are set to the last sample of the previous trial. For more detailed information, please refer to the vignette "Simulations of the Autocorrelated Bayesian Sampler".
    #' @param ... further arguments passed to the ABS model, see also `Details`.
    #' 
    #' @details
    #' The ABS model has two types of stopping rules: fixed and relative. The fixed stopping rule means that a fixed number of samples are drawn to complete the tasks such as estimations and confidence intervals. This rule applies to tasks such as estimation tasks. On the other hand, the relative stopping rule means that the model counts the difference in evidence between the two hypotheses, and terminates the sampling process whenever the accumulated difference exceeds a threshold. This rule applies to tasks such as two-alternative force choice tasks.
    #' 
    #' When the `stopping rule` is `"fixed"`, the following arguments are required:
    #' 
    #' - `n_sample` an integer of the fixed number of samples for each trial.
    #' - `trial_stim` a numeric vector of the stimulus of each trial.
    #' 
    #' When the `stopping rule` is `"relative"`, the following arguments are required:
    #' 
    #' - `delta` an integer of the relative difference between the number of samples supporting each hypothesis.
    #' - `dec_bdry` a numeric value of the decision boundary that separates the posterior hypothesis distribution.
    #' - `discrim` a numeric value of the stimuli discriminability.
    #' - `trial_stim` a factor that indicates the stimuli of each trial. It only consists of either one level or two levels. By definition, level 1 represents the stimulus below the decision boundary, while level 2 represents the stimulus above the decision boundary.
    #' - `prior_on_resp` a numeric vector for the Beta prior on responses. Defaults to `c(1,1)` representing the distribution `Beta(1,1)`.
    #' - `prior_depend` a boolean variable that control whether the prior on responses changes regarding the last stimulus. Defaults to `TRUE`. Please refer to the vignette for more information.
    #' - `max_iterations` an integer of the maximum length of the MC3 sampler. Defaults to 1000. The program will stop the sampling process after the length of the sampling sequence reaches to this limitation.
    #' 
    #' No values will be return after running this method, but the field `sim_results` will be updated instead. If the stopping rule is "fixed", `simulation_results` will be a data frame with five columns:
    #' \enumerate{
    #'  \item{trial: The index of trials;}
    #'  \item{samples: The samples of ABS sampler for the trial;}
    #'  \item{stimulus: The stimuli of the experiment;}
    #'  \item{rt: The response time;}
    #'  \item{point_est: The response of point estimation;}
    #'  }
    #' 
    #' On the other hand, if the stopping rule is "relative", `sim_results` will be a data frame with eight columns:
    #' \enumerate{
    #'  \item{trial: The index of trials;}
    #'  \item{samples: The samples of ABS sampler for the trial;}
    #'  \item{response: The response predicted by ABS;}
    #'  \item{stimulus: The stimuli of the experiment;}
    #'  \item{accuracy: Whether the response is the same as the feedback. 0 represents error, and 1 represents correct;}
    #'  \item{rt: The response time, including both the non-decision and the decision time;}
    #'  \item{confidence: The confidence of the response;}
    #'  \item{point_est: The response of point estimation.}
    #' }
    #' 
    #' @examples
    #' 
    #' trial_stim <- round(runif(5, 10, 50))
    #' zhuabs$simulate(stopping_rule='fixed', n_sample = 5, trial_stim = trial_stim)
    #' zhuabs$sim_results
    #' 
    #' zhuabs$reset_sim_results()
    #' trial_stim <- factor(sample(c('left', 'right'), 5, TRUE))
    #' zhuabs$simulate(stopping_rule='relative', 
    #'    delta = 4, dec_bdry = 0, 
    #'    discrim = 1, trial_stim = trial_stim
    #' )
    #' zhuabs$sim_results
    #' 
    simulate = function(stopping_rule, start_point=NA, ...){
      # Check samples
      if (is.data.frame(self$sim_results)){
        stop("Samples have been drawn. Please use the `reset_sim_results` method to reset the samples if you want to rerun the simulation.\n")
      }
      # Check stopping rule
      if (!(stopping_rule %in% c('fixed', 'relative'))){
        stop(paste0('The stopping rule "', stopping_rule, '" is not supported by ABS.'))
      }
      
      switch(
        stopping_rule,
        fixed = {
          private$stopping_rule <- 'fixed'
          private$simulate_fixed_sr(start_point = start_point, ...)
        },
        relative = {
          private$stopping_rule <- 'relative'
          private$simulate_relative_sr(start_point = start_point, ...)
        }
      ) # end of switch
      
      invisible(self)
    },
    
    #' @description
    #' This function calculates the confidence interval of the `simulate` method's results when the "fixed" stopping rule was used.
    #'
    #' @param conf_level the required confidence level.
    #'
    #' @details
    #' No values will be returned by this method. Instead, two new columns will be added to the `sim_results`:
    #' \enumerate{
    #'  \item{conf_interval_l: The lower bound of the confidence interval with the given level;}
    #'  \item{conf_interval_u: The upper bound of the confidence interval with the given level;}
    #'  }
    #'
    #' @examples
    #' zhuabs$confidence_interval(conf_level = 0.9)
    #'
    confidence_interval = function(conf_level){

      # Check sim_results
      if (!is.data.frame(self$sim_results)){
        stop("Please run the `estimate` method first.\n")
      }
      
      # Check the stopping rule
      if(private$stopping_rule == 'relative'){
        warning("The simulation results were geneated with the relative stopping rule.")
      }
      
      # Check conf_level  
      stopifnot('Argument "conf_level" should be a single value between 0 and 1.'=(is.numeric(conf_level) & length(conf_level) == 1 & conf_level >= 0 & conf_level <= 1))

      conf_interval <- t(sapply(self$sim_results$samples, function(samples) quantile(samples, probs = c((1-conf_level)/2, (1+conf_level)/2))))
      self$sim_results$conf_interval_l <- conf_interval[,1]
      self$sim_results$conf_interval_u <- conf_interval[,2]
      
      invisible(self)
    },
    
    
    #' @description
    #' This function is for resetting the `sim_results` to run new simulations.
    reset_sim_results = function(){
      self$sim_results <- NULL
      invisible(self)
    }
  ),
  
  private = list(
    
    simulate_fixed_sr = function(n_sample, trial_stim, start_point){
      
      #Check inputs
      stopifnot('Argument "n_sample" should be a single integer.'=is.numeric(n_sample))
      stopifnot('Argument "n_sample" should be a single integer.'=(n_sample %% 1==0 & length(n_sample) == 1))
      stopifnot('Argument "trial_stim" should be a numeric vector.'=(is.numeric(trial_stim)))
      
      # Check starting points
      if (any(!is.na(start_point))){
        if (any(is.na(start_point))) { # if start_point contains NA
          stop('Argument "start_point" contains NA.')
        }
        stopifnot('Argument "start_point" should be a numeric vector.' = (is.numeric(start_point)))
        stopifnot('The length of "start_point" should equal to the length of "trial_stim".' = (length(start_point) == length(trial_stim)))
      }
      
      # Prepare for the custom distributions
      n_trial <- length(trial_stim)
      if (is.null(self$custom_distr)){
        distr_name <- self$distr_name
        # place holders for custom distribution
        custom_distr <- replicate(n_trial, function(x){}, simplify = FALSE)
        custom_start <- 1
        use_custom <- FALSE
        
        # Check the inputs of `distr_params`
        if (length(self$distr_params) == 1){
          distr_add_params <- rep(self$distr_params, n_trial)
        } else {
          stopifnot('The length of "distr_params" should equal to either 1 or the length of "trial_stim".' = (length(self$distr_params) == n_trial))
          distr_add_params <- self$distr_params
        }
        
      } else {
        distr_name <- ""
        distr_add_params <- rep(1, n_trial)
        use_custom <- TRUE
        custom_distr <- self$custom_distr
        custom_start <- self$custom_start
      }
      
      # Start the simulation
      samples_fixed_sr <- Zhu23ABS_cpp(
        stop_rule_id = 1,
        trial_stim = trial_stim,
        distr_name = distr_name,
        distr_add_params = distr_add_params,
        custom_func = custom_distr,
        custom_start = custom_start,
        useCustom = use_custom,
        proposal_width = self$width,
        n_chains = self$n_chains,
        provided_start_point = start_point,
        stop_rule = n_sample,
        nd_time = self$nd_time,
        s_nd_time = self$s_nd_time,
        lambda = self$lambda
      )
      self$sim_results <- data.frame(do.call(rbind, samples_fixed_sr))
      rm(samples_fixed_sr)
    
      self$sim_results$point_est <- sapply(self$sim_results$samples, function(samples) samples[length(samples)])
      
      invisible(self)
    },
    
    simulate_relative_sr = function(delta, dec_bdry, discrim, trial_stim, start_point, prior_on_resp = c(1,1), prior_depend=TRUE, max_iterations=1000){
      
      #Check inputs
      stopifnot('Argument "delta" should be a single integer.'=is.numeric(delta))
      stopifnot('Argument "delta" should be a single integer.'=(delta %% 1==0 & length(delta) == 1))
      stopifnot('Argument "prior_on_resp" should be a numeric vector with two values.'=(is.numeric(prior_on_resp) && length(prior_on_resp) == 2))
      stopifnot('Argument "dec_bdry" should be a single numeric value.'=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
      stopifnot('Argument "discrim" should be a single numeric value.'=(is.numeric(discrim) && length(discrim) == 1))
      stopifnot('Argument "trial_stim" should be a factor.'=is.factor(trial_stim))
      stopifnot('Argument "prior_depend" should be a boolean variable.'=(isTRUE(prior_depend) || isFALSE(prior_depend)))
      stopifnot('Argument "max_iterations" should be a single numeric value.'=(is.numeric(max_iterations) && length(max_iterations) == 1))
      
      # Check starting points
      if (any(!is.na(start_point))){
        if (any(is.na(start_point))) { # if start_point contains NA
          stop('Argument "start_point" contains NA.')
        }
        stopifnot('Arguent "start_point" should be a numeric vector' = (is.numeric(start_point)))
        stopifnot('The length of "start_point" should equal to the length of "trial_stim".' = (length(start_point) == length(trial_stim)))
      }
      
      # Prepare for the custom distributions
      n_trial <- length(trial_stim)
      if (is.null(self$custom_distr)){
        distr_name <- self$distr_name
        # place holders for custom distributions
        custom_distr <- replicate(n_trial, function(x){}, simplify = FALSE)
        custom_start <- 1
        use_custom <- FALSE
        
        # Check the inputs of `distr_params`
        if (length(self$distr_params) == 1){
          distr_add_params <- rep(self$distr_params, n_trial)
        } else {
          stopifnot('The length of "distr_params" should equal to either 1 or the length of "trial_stim".' = (length(self$distr_params) == n_trial))
          distr_add_params <- self$distr_params
        }
        
      } else {
        distr_name <- ""
        distr_add_params <- rep(1, n_trial)
        use_custom <- TRUE
        custom_distr <- self$custom_distr
        custom_start <- self$custom_start
      }
      
      trial_stim_num <- as.numeric(trial_stim)
      stim_levels <- levels(trial_stim)
      stopifnot('Argument "trial_stim" should not have more than two levels.' = (length(stim_levels) <= 2))
      
      
      # Start the simulation
      samples_relative_sr <- Zhu23ABS_cpp(
        stop_rule_id = 2,
        trial_stim = trial_stim, 
        distr_name = distr_name,
        distr_add_params = distr_add_params,
        custom_func = custom_distr,
        custom_start = custom_start,
        useCustom = use_custom,
        proposal_width = self$width,
        n_chains = self$n_chains,
        provided_start_point = start_point,
        prior_on_resp = prior_on_resp,
        stop_rule = delta,
        nd_time = self$nd_time, 
        s_nd_time = self$s_nd_time,
        lambda = self$lambda,
        prior_depend = prior_depend,
        mc3_iterations = max_iterations,
        dec_bdry = dec_bdry, 
        discrim = discrim
      )
      self$sim_results <- data.frame(do.call(rbind, samples_relative_sr))
      rm(samples_relative_sr)
      
      self$sim_results$stimulus <- trial_stim
      self$sim_results$response <- stim_levels[as.numeric(self$sim_results$response)]
      self$sim_results$accuracy <- as.numeric(self$sim_results$accuracy)
      self$sim_results$rt <- as.numeric(self$sim_results$rt)
      self$sim_results$point_est <- sapply(self$sim_results$samples, function(samples) samples[length(samples)])
      self$sim_results$confidence <- as.numeric(self$sim_results$confidence)
      invisible(self)
    }
  )
  
)