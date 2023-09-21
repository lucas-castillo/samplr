#' ABS simulation
#'
#' This function is for simulating the Auto-correlated Bayesian Sampler.
#'
#' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
#' @param discrim The stimuli discriminability.
#' @param delta The stopping rule.
#' @param nd_time The non-decision time.
#' @param s_nd_time The range of the non-decision time.
#' @param lambda The rate parameter of the gamma distribution for decision time.
#' @param trial_fdbk The feedback of each trial. It should be either 0 (below the decision boundary) or 1 (above the decision boundary).
#' @param nChains The number of chains of the sampler.
#' @param width The proposal width of the sampler.
#' @param distr_name The type of the posterior hypothesis distribution.
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
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' simulation <- rABS(dec_bdry=0, discrim=1, delta=2, nd_time=0.3, s_nd_time=0.5, lambda=6, trial_fdbk=c(0, 0, 1, 0, 1), nChains=3, width=1)
#' 

rABS <- function(dec_bdry, discrim, delta, nd_time, s_nd_time, lambda, trial_fdbk, nChains, width, distr_name='norm') {
  # Check input
  start_point <- runif(nChains, min=-3, max=3) %>%
    as.matrix()
  mc3_iterations=100
  
  abs_sim <- ABS_sampler_cpp(
    start_point = start_point,
    trial_fdbk = trial_fdbk, 
    distr_name = distr_name, 
    mc3_iterations = mc3_iterations, 
    nChains = nChains, 
    dec_bdry = dec_bdry, 
    d_sepn = discrim, 
    delta = delta, 
    nd_time = nd_time, 
    s_nd_time = s_nd_time,
    er_lambda = lambda,
    proposal_width = width
  )
  
  x <- abs_sim %>% 
    do.call(cbind, .) %>% 
    t %>% 
    data.frame() %>% 
    tibble::tibble() %>% 
    tidyr::unnest(-c(samples, support))
  
  x$rept <- c(x$response[2:(length(x$response))] == x$response[1:(length(x$response) - 1)], NA)
  return(x)
}