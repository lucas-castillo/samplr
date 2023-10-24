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
//[[Rcpp::export]]
List sampler_mh_cpp(
    NumericVector start,
    NumericMatrix sigma_prop,
    int iterations,
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
  // Initialize variables ---------------------------------
  LogicalVector acceptances(iterations);
  int n_dim = start.size();
  dfunc pdf = managePDF(distr_name, distr_params, isMix, weights, false, custom_func, useCustom);
  
  NumericMatrix chain(iterations, n_dim);
  NumericMatrix proposals(iterations, n_dim);
  NumericMatrix jumps(iterations, n_dim);
  NumericMatrix true_jumps(iterations, n_dim);
  NumericMatrix ps(1, iterations);
  
  // first row is start
  chain.row(0) = start;
  ps(0,0) = pdf(start);
  
  
  // Run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){
    NumericVector accept;
    if (i == 1){
      accept = autocorrelated_metropolis_step_cpp(
        chain,      // NumericMatrix &chain, 
        proposals, // NumericMatrix &proposals, 
        jumps, // NumericMatrix &jumps, 
        true_jumps, // NumericMatrix &true_jumps, 
        i,        // const int &currentIndex, 
        ps(0,i-1), // const double &lastP, 
        sigma_prop, // const NumericMatrix &sigma_prop, 
        pdf, // dfunc &pdf, 
        discreteValues, // const bool &discreteValues, 
        1, // const double &beta, 
        0 // const double &alpha
      );
    } else{
      accept = autocorrelated_metropolis_step_cpp(
        chain,      // NumericMatrix &chain, 
        proposals, // NumericMatrix &proposals, 
        jumps, // NumericMatrix &jumps, 
        true_jumps, // NumericMatrix &true_jumps, 
        i,        // const int &currentIndex, 
        ps(0,i-1), // const double &lastP, 
        sigma_prop, // const NumericMatrix &sigma_prop, 
        pdf, // dfunc &pdf, 
        discreteValues, // const bool &discreteValues, 
        1, // const double &beta, 
        alpha // const double &alpha
      );
    }
    
    ps(0,i) = accept(0);
    acceptances(i) = (bool)(accept(1));
  }
  
  return List::create(chain, proposals, acceptances, ps, jumps, true_jumps);
}
