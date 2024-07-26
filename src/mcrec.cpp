// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


#include "hmc_utils.h"
#include "pdf_manage.h"
using namespace Rcpp;

///'@export
// [[Rcpp::export]]
List sampler_mc_rec_cpp(
    NumericVector start,
    int nChains, // for non-tempered versions of HMC and HOR, set this to one
    double delta_T, // for all chains having the same temp, set to 0
    bool swap_all,
    double iterations,
    StringVector distr_name,
    List distr_params,
    bool discreteValues,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom,
    double epsilon,
    int L,
    double alpha // for HMC, set this to zero
)
{
  if (alpha < -1 || alpha > 1){
    stop("Alpha must be between -1 and +1 (inclusive)");
  }
  
  // init vars
  NumericVector acceptances(nChains); // counts accept per chain
  int swap_attempts = 0; // counter
  int swap_accepts= 0; // counter
  NumericMatrix swaps(iterations*nChains, 3); //swap History. Columns: iteration, chain1, chain2
  int n_dim = start.length();
  
  NumericMatrix chain(iterations * nChains, n_dim);
  NumericMatrix proposals(iterations*nChains, n_dim);
  NumericVector beta(nChains); // vector of temperatures
  NumericMatrix ps(nChains, iterations); // densities
  NumericMatrix momentumsHistory(iterations*nChains, n_dim);
  NumericMatrix momentumsUse(iterations*nChains, n_dim);
  NumericVector momentum;
  
  dfunc log_pdf = managePDF(distr_name, distr_params, isMix, weights, true, custom_func, useCustom);
  // init the chain information
  for (int i = 0; i < nChains; i++){
    chain.row(0 + iterations * i) = start;
    beta(i) = (1 + delta_T * i);
  }
  
  int nSwaps; // how many swaps per iteration?
  if (swap_all){
    nSwaps = floor(nChains/2);
  } else {
    nSwaps = 1;
  }
  
  // v holds the chain ids so that they are swapped
  NumericVector v(nChains);
  v = seq(0, nChains - 1);
  
  // run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){
    
    if (i % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    for (int ch = 0; ch < nChains; ch++){
      if (i == 1){
        momentum = drawMomentum(n_dim) * beta(ch);
      } else{
        momentum = momentumsUse.row((i + iterations * ch) - 1);
        momentum = RecycledMomentumUpdate(momentum, alpha, beta(ch));
      }
      
      NumericVector theta_prime  = chain.row((i + iterations * ch) - 1);
      NumericVector momentum_prime = clone(momentum);
      
      // leapfrog for each L step
      leapfrog_step_cpp(theta_prime, momentum_prime, epsilon, log_pdf, L, beta(ch));
      
      // Metropolis - Hastings Acceptance, using the joint density of position + momentum
      
      double top =  exp(joint_d(theta_prime, momentum_prime, log_pdf, beta(ch)));
      double bottom =  exp(joint_d(chain.row((i + iterations * ch) - 1), momentum, log_pdf, beta(ch)));
      
      double alpha = top/bottom;
      if (R::runif(0,1) <= alpha){
        chain.row(i + iterations * ch) = theta_prime;
        momentumsUse.row(i + iterations * ch) = -1 * momentum_prime;
        momentumsHistory.row(i + iterations * ch) = -1 * momentum_prime;
        acceptances(ch) = acceptances(ch) + 1;
      } else{
        chain.row(i + iterations * ch) = chain.row((i + iterations * ch) -1);
        momentumsUse.row(i + iterations * ch) = -1 * momentumsUse.row((i + iterations * ch) -1);
        momentumsHistory.row(i + iterations * ch) = -1 * momentumsUse.row((i + iterations * ch) -1);
      }
    }
    
    // once a step in every chain has been done, proceed to swap chains
    if (nChains > 1){
      // arrange chains randomly
      
      v = sample(v, nChains, false);
      
      // swap nSwaps times (depending on swap_all)
      for (int k = 0; k < nSwaps; k++){
        swap_attempts++;
        int m = v[k*2];
        int n = v[k*2 + 1];
        // chains are swapped with probability alpha, which is the ratio between:
        // - the product of the density of each chain's location at the temperature of the other chain, and
        // - the product of the density of each chain's location at their own temperatures
        
        double m_pdf = exp(joint_d(chain.row(i + iterations * m), momentumsUse.row(i + iterations * m), log_pdf, beta(m)));
        double n_pdf = exp(joint_d(chain.row(i + iterations * n), momentumsUse.row(i + iterations * n), log_pdf, beta(n)));
        
        double m_pdf_swapped = exp(joint_d(chain.row(i + iterations * m), momentumsUse.row(i + iterations * m), log_pdf, beta(n)));
        double n_pdf_swapped = exp(joint_d(chain.row(i + iterations * n), momentumsUse.row(i + iterations * n), log_pdf, beta(m)));        // double n_pdf_swapped = exp(joint_d(chain.row(i + iterations * n), momentums.row(i + iterations * n), log_pdf, beta(m)));
        
        double top = m_pdf_swapped * n_pdf_swapped; // pow(m_pdf, beta(n)) * pow(n_pdf, beta(m));
        double bottom = m_pdf * n_pdf; //pow(m_pdf,beta(m)) * pow(n_pdf, beta(n));
        
        
        if ((bottom != 0 && R::runif(0,1) <= top/bottom) || (bottom == 0 && top > 0)){
          
          NumericVector temp = chain.row(i + iterations * m);
          chain.row(i + iterations * m) = chain.row(i + iterations * n);
          chain.row(i + iterations * n) = temp;
          
          NumericVector tempMomentum_m = momentumsUse.row(i + iterations * m);
          NumericVector tempMomentum_n = momentumsUse.row(i + iterations * n);
          
          if (i != (iterations-1)){
            momentumsUse.row(i + iterations * m) = drawMomentum(n_dim) * beta(m);
            momentumsUse.row(i + iterations * n) = drawMomentum(n_dim) * beta(n);
          }
          
          swaps.row(swap_accepts) = NumericVector::create(i+1, m, n);
          swap_accepts++;
        }
      }
    }
  }
  
  
  return List::create(
    _["chain"] = chain,
    _["proposals"] = proposals,
    _["beta"] = beta,
    _["swaps"] = swaps,
    _["swap_accepts_attempts"] = NumericVector::create(swap_accepts, swap_attempts),
    _["acceptances"] = acceptances,
    _["momentums"] = momentumsUse,
    _["momentumsHistory"] = momentumsHistory
  );
}
