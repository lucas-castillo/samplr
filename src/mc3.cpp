// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

#include "mh_utils.h"
#include "pdf_manage.h"
using namespace Rcpp;

///'@export
// [[Rcpp::export]]
List sampler_mc3_cpp(
    NumericMatrix start, // Numeric Matrix of starts for each chain.
    int nChains,
    NumericMatrix sigma_prop,
    double delta_T,
    bool swap_all,
    double iterations,
    StringVector distr_name,
    List distr_params,
    bool discreteValues,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom,
    double alpha=0
)
{
  NumericVector acceptances(nChains);
  int swap_attempts = 0;
  int swap_accepts= 0;
  NumericMatrix swaps(iterations*nChains, 3);
  int n_dim = start.ncol();
  
  NumericMatrix chain(iterations*nChains, n_dim);
  NumericMatrix proposals(iterations*nChains, n_dim);
  NumericVector beta(nChains);
  NumericMatrix ps(nChains, iterations);
  NumericMatrix jumpsHistory(iterations*nChains, n_dim);
  NumericMatrix jumpsUse(iterations*nChains, n_dim);
  NumericVector momentum;
  
  dfunc pdf = managePDF(distr_name, distr_params, isMix, weights, false, custom_func, useCustom);
  for (int i = 0; i < nChains; i++){
    chain.row(0 + iterations * i) = start.row(i);
    
    
    beta(i) = 1 / (1 + delta_T * i);
    ps(i,0) = pdf(start.row(i));
  }
  int nSwaps;
  if (swap_all){
    nSwaps = floor(nChains/2);
  } else {
    nSwaps = 1;
  }
  
  NumericVector v(nChains);
  v = seq(0, nChains - 1);
  
  // run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){

    for (int ch = 0; ch < nChains; ch++){
      NumericVector accept;
      if (i ==1){
        accept =  autocorrelated_metropolis_step_cpp(
          chain, // NumericMatrix &chain, 
          proposals, // NumericMatrix &proposals, 
          jumpsUse, // NumericMatrix &jumps, 
          jumpsHistory,// NumericMatrix &true_jumps, 
          i + iterations * ch, // const int &currentIndex, 
          ps(ch, i-1), // const double &lastP, 
          sigma_prop, // const NumericMatrix &sigma_prop, 
          pdf, // dfunc &pdf, 
          discreteValues, // const bool &discreteValues, 
          beta(ch), // const double &beta, 
          0  // const double &alpha
        );
      } else{
        accept =  autocorrelated_metropolis_step_cpp(
          chain, // NumericMatrix &chain, 
          proposals, // NumericMatrix &proposals, 
          jumpsUse, // NumericMatrix &jumps, 
          jumpsHistory,// NumericMatrix &true_jumps, 
          i + iterations * ch, // const int &currentIndex, 
          ps(ch, i-1), // const double &lastP, 
          sigma_prop, // const NumericMatrix &sigma_prop, 
          pdf, // dfunc &pdf, 
          discreteValues, // const bool &discreteValues, 
          beta(ch), // const double &beta, 
          alpha // const double &alpha
        );
      }
      
      ps(ch, i) = accept(0);
      acceptances(ch) = acceptances(ch) + accept(1);
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
        NumericVector t1 = chain.row(i+iterations * m);
        NumericVector t2 = chain.row(i+iterations * n);
        double m_pdf = pdf(chain.row(i + iterations * m));
        double n_pdf = pdf(chain.row(i + iterations * n));
        
        double top = pow(m_pdf, beta(n)) * pow(n_pdf, beta(m));
        double bottom = pow(m_pdf,beta(m)) * pow(n_pdf, beta(n));
        
        
        if ((bottom != 0 && R::runif(0,1) <= top/bottom) || (bottom == 0 && top > 0)){
          
          // Swap Positions
          NumericVector temp = chain.row(i + iterations * m);
          chain.row(i + iterations * m) = chain.row(i + iterations * n);
          chain.row(i + iterations * n) = temp;
          // Record Swap
          swaps.row(swap_accepts) = NumericVector::create(i+1, m, n);
          swap_accepts++;
          // Swap probabilities
          t1 = chain.row(i+iterations * m);
          t2 = chain.row(i+iterations * n);
          double TEMP = ps(m, i);
          ps(m,i) = ps(n, i);
          ps(n,i) = TEMP;
          // Swap last jump
          if (i != (iterations-1)){
            // Reset jump after swap
            NumericVector zeros (n_dim);
            arma::mat random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
            NumericVector random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
            
            jumpsUse.row(i + iterations * m) = random_jump;
            
            
            random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
            random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
            
            jumpsUse.row(i + iterations * n) = random_jump;
          }
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
    _["jumps"] = jumpsUse,
    _["jumpsHistory"] = jumpsHistory
  );
}

