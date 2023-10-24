// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// 
#include "pdf_manage.h"


using namespace Rcpp;

NumericVector alpha_trick(
    NumericVector random_jump,
    NumericVector last_jump,
    double alpha
){
  return alpha * last_jump + pow((1 - pow(alpha, 2)), .5) * random_jump;
}



NumericVector autocorrelated_metropolis_step_cpp(
    NumericMatrix &chain,
    NumericMatrix &proposals,
    NumericMatrix &jumps,
    NumericMatrix &true_jumps,
    const int &currentIndex,
    const double &last_prob,
    const NumericMatrix &sigma_prop,
    dfunc &pdf,
    const bool &discreteValues,
    const double &beta,
    const double &alpha
){
  NumericVector current_x = chain.row(currentIndex - 1);
  NumericVector last_jump = jumps.row(currentIndex - 1);
  NumericVector zeros (current_x.size());
  arma::mat random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
  NumericVector random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
  
  NumericVector next_jump = alpha_trick(random_jump, last_jump, alpha);
  true_jumps.row(currentIndex) = next_jump;
  NumericVector proposal = current_x + next_jump;
  
  // if dist is discrete round the proposal to nearest int
  if (discreteValues){
    for (int i = 0; i < proposal.length(); i++){
      proposal(i) = round(proposal(i));
    }
  }
  
  // update proposals matrix
  proposals.row(currentIndex) = proposal;
  
  // calculate proposal density
  double prob_prop = pdf(proposal);
  
  // proposal is accepted with probability prob_prop / prob_curr
  if (last_prob != 0){
    double ratio = prob_prop / last_prob;
    
    // The beta parameter (temperature), beta <= 1,
    // increases the value of the ratio making hotter chains
    // more likely to accept proposals
    if ((ratio >= 1) || (R::runif(0,1) < pow(ratio, beta))){
      chain.row(currentIndex) = proposal;
      jumps.row(currentIndex) = next_jump;
      return NumericVector::create(prob_prop, 1);
      
    }
  } else if (prob_prop > 0) {
    chain.row(currentIndex) = proposal;
    jumps.row(currentIndex) = next_jump;
    return NumericVector::create(prob_prop, 1);
  }
  chain.row(currentIndex) = current_x;
  jumps.row(currentIndex) = next_jump * -1;
  return NumericVector::create(last_prob, 0);
}
