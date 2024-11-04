#' Z Identities
#' 
#' Calculates identities Z1 to Z18 as defined in \insertCite{costello2016PeopleConditionalProbability,zhu2020BayesianSamplerGeneric}{samplr}. Probability theory predicts that these will all equal 0. 
#' 
#' If some of the probability estimates are not given, calculation will proceed and equalities that cannot be calculated will be coded as NA. 
#' 
#'
#' @param a,b,a_and_b,a_or_b,a_given_b,b_given_a,a_given_not_b,b_given_not_a,a_and_not_b,b_and_not_a Probability estimates given by participants
#' @param not_a,not_b Probability estimates given by participants. If not given, they'll default to 1-a and 1-b respectively
#'
#' @return Dataframe with identities Z1 to Z18 
#' @export
#'
#' @examples
#'Z_identities(
#'  a=.5, 
#'  b=.1, 
#'  a_and_b=.05, 
#'  a_or_b=.55, 
#'  a_given_b=.5,
#'  b_given_a=.1,
#'  a_given_not_b=.5,
#'  b_given_not_a=.1,
#'  a_and_not_b=.45,
#'  b_and_not_a=.05,
#'  )
#'#Get identities for a set of participants
#'library(magrittr)
#'library(dplyr)
#'library(tidyr)
#'data.frame(
#'  ID = LETTERS[1:20],
#'  a=runif(20),
#'  b=runif(20),
#'  a_and_b=runif(20),
#'  a_or_b=runif(20),
#'  a_given_b=runif(20),
#'  b_given_a=runif(20),
#'  a_given_not_b=runif(20),
#'  b_given_not_a=runif(20),
#'  a_and_not_b=runif(20),
#'  b_and_not_a=runif(20),
#'  not_a=runif(20),
#'  not_b=runif(20)
#') %>% 
#'  group_by(ID) %>% 
#'  do(
#'    Z_identities(
#'      .$a,
#'      .$b,
#'      .$a_and_b,
#'      .$a_or_b,
#'      .$a_given_b,
#'      .$b_given_a,
#'      .$a_given_not_b,
#'      .$b_given_not_a,
#'      .$a_and_not_b,
#'      .$b_and_not_a,
#'      .$not_a,
#'      .$not_b
#'    )
#'  )
Z_identities <- function(
    a=NULL, 
    b=NULL, 
    a_and_b=NULL, 
    a_or_b=NULL, 
    a_given_b=NULL,
    b_given_a=NULL,
    a_given_not_b=NULL,
    b_given_not_a=NULL,
    a_and_not_b=NULL,
    b_and_not_a=NULL,
    not_a=NULL,
    not_b=NULL
    ){
  # Check inputs
  is_prob <- function(x){
    x <- x[!is.na(x)]
    if (any(x<0)) stop("Probabilites cannot be negative")
    if (any(x>1)) stop("Probabilites cannot be larger than 1")
  }
  for (p in c(a,
             b,
             not_a,
             not_b,
             a_and_b,
             a_or_b,
             a_given_b,
             b_given_a,
             a_given_not_b,
             b_given_not_a,
             a_and_not_b,
             b_and_not_a)){
    if (!is.null(p)){is_prob(p)}
  }
  
  # add complements if not given
  if (is.null(not_a)){
    not_a = 1-a
  }
  if (is.null(not_b)){
    not_b = 1-b
  }
  
  safe_eval <- function(exp){
    if (length(eval(exp)) == 0) return(NA) else return(round(eval(exp),10))
  }
  
  data.frame(
    z1 = safe_eval(expression(a + b - a_and_b - a_or_b)),
    z2 = safe_eval(expression(a + b_and_not_a - b - a_and_not_b)),
    z3 = safe_eval(expression(a + b_and_not_a - a_or_b)),
    z4 = safe_eval(expression(b + a_and_not_b - a_or_b)),
    z5 = safe_eval(expression(a_and_not_b + a_and_b - a)),
    z6 = safe_eval(expression(b_and_not_a + a_and_b - b)),
    z7 = safe_eval(expression(a_and_not_b + b_and_not_a + a_and_b - a_or_b)),
    z8 = safe_eval(expression(a_and_not_b + b_and_not_a + 2*a_and_b - a - b)),
    z9 = safe_eval(expression(a_given_b*b - b_given_a*a)),
    z10= safe_eval(expression(a_given_b*b + a_given_not_b*not_b - a)),
    z11= safe_eval(expression(b_given_a*a + b_given_not_a*not_a - b)),
    z12= safe_eval(expression(b_given_a*a + a_given_not_b*not_b - a)),
    z13= safe_eval(expression(a_given_b*b + b_given_not_a*not_a - b)),
    z14= safe_eval(expression(a_given_not_b*not_b + b - b_given_not_a*not_a - a)),
    z15= safe_eval(expression(a_and_b - a_given_b*b)),
    z16= safe_eval(expression(a_and_b - b_given_a*a)),
    z17= safe_eval(expression(a_and_b - a + a_given_not_b*not_b)),
    z18= safe_eval(expression(a_and_b - b + b_given_not_a*not_a))
  )
}

get_true_probabilities <- function(
    a_and_b,
    b_and_not_a,
    a_and_not_b,
    not_a_and_not_b
  ){
  # +-------+-------------+-----------------+
  # |   .   |      a      |      not_a      |
  # +-------+-------------+-----------------+
  # | b     | a_and_b     | b_and_not_a     |
  # | not_b | a_and_not_b | not_a_and_not_b |
  # +-------+-------------+-----------------+
  # Normalize probabilities
  base = a_and_b + b_and_not_a + a_and_not_b + not_a_and_not_b
  a_and_b = a_and_b / base
  b_and_not_a = b_and_not_a / base
  a_and_not_b = a_and_not_b / base
  not_a_and_not_b = not_a_and_not_b / base
  # Get the other probabilities
  a = a_and_b + a_and_not_b
  b = a_and_b + b_and_not_a 
  not_a = 1 - a
  not_b = 1 - b
  a_or_b = 1 - not_a_and_not_b
  a_or_not_b = 1 - a_and_not_b
  b_or_not_a = 1- b_and_not_a
  not_a_or_not_b = 1 - a_and_b
  # Conditional Probabilities
  a_given_b = a_and_b / (a_and_b + b_and_not_a)
  not_a_given_b = b_and_not_a / (a_and_b + b_and_not_a)
  a_given_not_b = a_and_not_b / (a_and_not_b + not_a_and_not_b)
  not_a_given_not_b = not_a_and_not_b / (a_and_not_b + not_a_and_not_b)
  b_given_a = a_and_b / (a_and_b + a_and_not_b)
  not_b_given_a = a_and_not_b / (a_and_b + a_and_not_b)
  b_given_not_a = b_and_not_a / (b_and_not_a + not_a_and_not_b)
  not_b_given_not_a = not_a_and_not_b / (b_and_not_a + not_a_and_not_b)
  
  return(
    list(
      a_and_b = a_and_b,
      b_and_not_a = b_and_not_a,
      a_and_not_b = a_and_not_b,
      not_a_and_not_b = not_a_and_not_b,
      a = a,
      b = b,
      not_a = not_a,
      not_b = not_b,
      a_or_b = a_or_b,
      a_or_not_b = a_or_not_b,
      b_or_not_a = b_or_not_a,
      not_a_or_not_b = not_a_or_not_b,
      a_given_b = a_given_b,
      not_a_given_b = not_a_given_b,
      a_given_not_b = a_given_not_b,
      not_a_given_not_b = not_a_given_not_b,
      b_given_a = b_given_a,
      not_b_given_a = not_b_given_a,
      b_given_not_a = b_given_not_a,
      not_b_given_not_a = not_b_given_not_a
    )
  )
}

#' Bayesian Sampler Model 
#' 
#' As described in \insertCite{zhu2020BayesianSamplerGeneric}{samplr}. Vectors can be provided for each parameter, allowing multiple estimates at once. 
#'
#' @param a_and_b,b_and_not_a,a_and_not_b,not_a_and_not_b True probabilites for the conjuctions and disjunctions of A and B. Must add to 1.
#' @param beta Prior parameter.
#' @param N Number of samples drawn
#' @param N2 Optional. Number of samples drawn for conjunctions and disjunctions. (called N' in the paper). If not given, it will default to N2=N. Must be equal or smaller than N. 
#' @param return Optional. Either "mean", "variance" or "simulation". 

#' @return If return="mean" or return="variance", named list with predicted probabilities for every possible combination of A and B, or the expected variance of those predictions. If return="simulation", simulated predictions instead. 
#' @export
#'
#' @examples
#' Bayesian_Sampler(
#'     a_and_b = c(.4, .25),
#'     b_and_not_a = c(.4,  .25),
#'     a_and_not_b = c(.1, .25),
#'     not_a_and_not_b = c(.1, .25),
#'     beta = 1,
#'     N <- c(10, 12),
#'     N2 <- c(10, 10)
#' )

Bayesian_Sampler <- function(
    a_and_b,
    b_and_not_a,
    a_and_not_b,
    not_a_and_not_b,
    beta, N, N2=NULL, 
    return="mean"){
  
  if (sd(
    lengths(
      list(a_and_b,b_and_not_a,a_and_not_b,not_a_and_not_b)
    )
  ) != 0){stop("Probability vectors must all be the same length")}
  if (length(beta) != 1 & length(beta) != length(a_and_b)){
    stop("Beta length should be either 1 or the length of the probability vector.")
  }
  if (length(N) != 1 & length(N) != length(a_and_b)){
    stop("N length should be either 1 or the length of the probability vector.")
  }
  if (is.null(N2)){N2 <- N} else if (length(N2) != 1 & length(N2) != length(a_and_b)){
    stop("N2 length should be either 1 or the length of the probability vector.")
  } else if (any(N2 > N)){
    stop("N2 is larger than N. N2 <= N is required.")
  }
  
  sums <- apply(matrix(c(a_and_b,
                         b_and_not_a,
                         a_and_not_b,
                         not_a_and_not_b), ncol = 4), 1, sum)
  
  if(!all.equal(sums, rep(1, length(sums)))){stop("Probabilities must add up to 1")}
  
  true_probabilities <- get_true_probabilities(
    a_and_b, 
    b_and_not_a, 
    a_and_not_b, 
    not_a_and_not_b
  )
  get_mean <- function(p, N, beta){
    p*N/(N+2*beta)+beta/(N+2*beta)
  }
  predicted_means <- list()
  for (name in names(true_probabilities)){
      if (name %in% c( # treat conjunctions differently
        "b_or_not_a", "not_a_or_not_b", "a_or_b", "a_or_not_b", 
        "a_and_b", "b_and_not_a", "a_and_not_b", "not_a_and_not_b")){
        predicted_means[[name]] <- get_mean(true_probabilities[[name]], N2, beta)
      } else{
        predicted_means[[name]] <- get_mean(true_probabilities[[name]], N, beta)
      }
  }
  return(predicted_means)
}

#' Mean Variance Estimates
#' 
#' Estimates number of samples and prior parameters of the Bayesian Sampler using the Mean/Variance relationship as shown by \insertCite{sundh2023UnifiedExplanationVariability}{samplr}. For consistency with the Bayesian Sampler function we call beta the prior parameter, and b0 and b1 slope and intercept respectively. 
#'
#' @param rawData Dataframe with the following column variables for N repetitions of each unique query: participant ID ('id'), response query 1, response query 2, ... , response query N
#' @param idCol Name of the 'ID' column.
#'
#' @return A dataframe with values for the intercept (b0) and slope (b1) of the estimated regression, as well as estimates for N, d, and beta (termed b in the paper) for each participant. 
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(magrittr)
#' library(samplrData)
#' data <- sundh2023.meanvariance.e3 %>%
#'   group_by(ID, querydetail) %>% 
#'   mutate(iteration = LETTERS[1:n()]) %>% 
#'   pivot_wider(id_cols = c(ID, querydetail), 
#'       values_from = estimate, names_from = iteration) %>% 
#'   mutate(across(where(is.numeric), \(x){x/100})) %>% 
#'   ungroup %>% 
#'   select(-querydetail)
#' head(data)
#' head(Mean_Variance(data, "ID"))
Mean_Variance <- function(rawData, idCol){
  #Step 1: Prepare data
  rawDataMatrix <- data.matrix(rawData[colnames(rawData)[colnames(rawData) != idCol]])
  meansByQuery <- apply(rawDataMatrix, 1, mean)
  varianceByQuery <- apply(rawDataMatrix, 1, var)
  preparedData <- cbind(rawData[idCol], varianceByQuery, meansByQuery*(1-meansByQuery))
  colnames(preparedData) <- c('id','varia','expect')
  #Step 2: Apply frequentist regression model
  linModel <- lme4::lmer(varia ~ expect + (expect | id), data = preparedData)
  
  coefficients <- coef(linModel)$id
  colnames(coefficients) <- c("b0", "b1")
  
  coefficients$N <- 1 / coefficients$b1
  coefficients$d <- (1 - sqrt(coefficients$N * coefficients$b0 * 4 + 1)) / 2
  coefficients$beta <- (coefficients$N * coefficients$d)  / (1 - 2 * coefficients$d)
  rownames(coefficients) <- NULL

  return(cbind(rawData[idCol], coefficients))                   
}
  
