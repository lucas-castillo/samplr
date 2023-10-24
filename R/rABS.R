#' ABS simulation
#'
#' This function is for simulating the Auto-correlated Bayesian Sampler \insertCite{zhuAutocorrelatedBayesian2023}{samplr}.
#'
#' @param dec_bdry The decision boundary that separates the posterior hypothesis distribution.
#' @param discrim The stimuli discriminability.
#' @param delta The threshold of the relative difference. It should be an integer.
#' @param nd_time The non-decision time.
#' @param s_nd_time The range of the non-decision time.
#' @param lambda The rate parameter of the gamma distribution for decision time.
#' @param trial_fdbk The feedback of each trial. It should be a character only consist of 0 (below the decision boundary) and/or 1 (above the decision boundary).
#' @param nChains The number of chains of the sampler. It should be an integer.
#' @param width The proposal width of the sampler.
#' @param distr_name The type of the posterior hypothesis distribution.
#' @param mc3_iterations The number of iterations for each chunk. See details for more information.
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
#' @details
#' `rABS` function runs by chunks. If the sampling sequence reach the stopping rule inside a chunk, then the sampler of will cut the sequence at that point and start a new trial simulation. Otherwise, it will begin a new chunk at the end of the previous one and reiterate the whole process.
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
#' simulation <- rABS(dec_bdry=0, discrim=1, delta=2, nd_time=0.3, s_nd_time=0.5, lambda=6, trial_fdbk=c(0, 0, 1, 0, 1), nChains=3, width=1)
#' 

rABS <- function(dec_bdry, discrim, delta, nd_time, s_nd_time, lambda, trial_fdbk, nChains, width, distr_name='norm', mc3_iterations=100) {
  
  # Check missing input
  if(any(missing(dec_bdry), missing(discrim), missing(delta), missing(nd_time), missing(s_nd_time),
         missing(lambda), missing(trial_fdbk), missing(nChains), missing(width))) stop("dec_bdry, discrim, delta, nd_time, s_nd_time, lambda, trial_fdbk, nChains and/or width must be supplied")
  # Check variable types
  stopifnot("dec_bdry should be a single numeric value."=(is.numeric(dec_bdry) && length(dec_bdry) == 1))
  stopifnot("discrim should be a single numeric value."=(is.numeric(discrim) && length(discrim) == 1))
  stopifnot("delta should be an integer."=(delta%%1 == 0))
  stopifnot("nd_time should be a single numeric value."=(is.numeric(nd_time) && length(nd_time) == 1))
  stopifnot("s_nd_time should be a single numeric value."=(is.numeric(s_nd_time) && length(s_nd_time) == 1))
  stopifnot("lambda should be a single numeric value."=(is.numeric(lambda) && length(lambda) == 1))
  stopifnot("trial_fdbk should be a vector with either 0 and/or 1."=all(trial_fdbk %in% c(0, 1)))
  stopifnot("nChains should be an integer."=(nChains%%1 == 0))
  stopifnot("width should be a single numeric value."=(is.numeric(width) && length(width) == 1))

  
  start_point <- runif(nChains, min=-3, max=3) %>%
    as.matrix()
  
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
  
  if (length(trial_fdbk)>1){
    x$rept <- c(NA, x$response[2:(length(x$response))] == x$response[1:(length(x$response) - 1)])
  } else {
    x$rept <- NA
  }
  return(x)
}