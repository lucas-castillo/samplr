#' ABS simulation
#' 
#' @description
#' This function is for simulating the Auto-correlated Bayesian Sampler \insertCite{@ABS, @zhuAutocorrelatedBayesian2023}{samplr}. It is a sequential sampling model assuming people draw autocorrelated samples from memory or beliefs, i.e., posterior of hypotheses.
#' 
#' @details
#' The `rABS` function in R runs simulations in chunks of fixed length determined by the `mc3_iterations` argument. In each trial simulation, the function generates the first chunk and checks if the simulated sequence meets the stopping rule within that chunk. If the stopping rule is met, the function cuts the chunk at that point, returns the sequence, and starts a new trial simulation. If the stopping rule is not met, the function proceeds to the next chunk, beginning where the previous chunk ended, and repeats the entire process.
#' 
#'
#' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
#' @param discrim The stimuli discriminability.
#' @param delta The threshold of the relative difference. It should be an integer.
#' @param prior_on_resp The beta prior on responses
#' @param nd_time The non-decision time.
#' @param s_nd_time The range of the non-decision time.
#' @param lambda The rate parameter of the gamma distribution for decision time.
#' @param trial_fdbk The feedback of each trial. It should be a numeric vector only consist of 0 (below the decision boundary) and/or 1 (above the decision boundary).
#' @param n_chains The number of chains of the sampler. It should be an integer.
#' @param width The proposal width of the sampler.
#' @param distr_name The type of the posterior hypothesis distribution.
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
#' 
#' @importFrom magrittr %>%
#' @importFrom Rdpack reprompt
#' 
#' @references
#'    \insertRef{zhuAutocorrelatedBayesian2023}{samplr}
#' 
#' @export
#'
#' @examples
#' simulation <- rABS(dec_bdry=0, discrim=1, delta=5, prior_on_resp = c(1, 0), nd_time=0.3, s_nd_time=0.5, lambda=6, trial_stim=c(0, 0, 1, 0, 1), n_chains=3, width=1)
#' 

rABS <- function(dec_bdry, discrim, delta, prior_on_resp, nd_time, s_nd_time, lambda, trial_stim, n_chains, width, distr_name='norm', mc3_iterations=100) {
  
  # Check variable types
  stopifnot("dec_bdry should be a single numeric value."=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
  stopifnot("discrim should be a single numeric value."=(is.numeric(discrim) && length(discrim) == 1))
  stopifnot("delta should be an integer."=(delta%%1 == 0))
  stopifnot("prior_on_resp should be a numeric vector with two values."=(is.numeric(prior_on_resp) && length(prior_on_resp) == 2))
  stopifnot("nd_time should be a single numeric value."=(is.numeric(nd_time) && length(nd_time) == 1))
  stopifnot("s_nd_time should be a single numeric value."=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
  stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
  stopifnot("trial_stim should be a vector with either 0 and/or 1."=all(trial_stim %in% c(0, 1)))
  stopifnot("n_chains should be an integer."=(n_chains%%1 == 0))
  stopifnot("width should be a single numeric value."=(is.numeric(width) && length(width) == 1))
  stopifnot("mc3_iterations should be be an integer"=(mc3_iterations%%1 == 0))
  
  start_point <- stats::runif(n_chains, min=-3, max=3) %>%
    as.matrix()
  
  tafc_sim <- Zhu23ABS_cpp(
    task_id = 2,
    start_point = start_point,
    trial_stim = trial_stim, 
    distr_name = distr_name, 
    mc3_iterations = mc3_iterations, 
    n_chains = n_chains, 
    dec_bdry = dec_bdry, 
    discrim = discrim,
    prior_on_resp = prior_on_resp,
    stop_rule = delta, 
    nd_time = nd_time, 
    s_nd_time = s_nd_time,
    er_lambda = lambda,
    proposal_width = width
  )
  tafc_df <- data.frame(do.call(rbind, tafc_sim))
  return(tafc_df)
}