// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

#include "hmc_utils.h"
#include "pdf_manage.h"
using namespace Rcpp;

///'@export
// [[Rcpp::export]]
List sampler_hmc_cpp(
    NumericVector start,
    StringVector distr_name,
    List distr_params,
    double epsilon,
    int L,
    int iterations,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom
)
{
  // init vars
  dfunc log_pdf = managePDF(distr_name, distr_params, isMix, weights, true, custom_func, useCustom);
  int dim = start.size();
  NumericMatrix chain(iterations, dim);
  NumericMatrix momentums(iterations, dim);
  int acceptances = 0;
  chain.row(0) = start;
  NumericVector momentum;
  
  for (int i = 1; i < iterations; i++){
    // draw a sample of momentum
    
    momentum = drawMomentum(dim);
    
    // initialize vars
    NumericVector theta_prime  = chain.row(i-1);
    NumericVector momentum_prime = clone(momentum);
    
    // leapfrog for each L step
    leapfrog_step_cpp(theta_prime, momentum_prime, epsilon, log_pdf, L);
    
    // Metropolis - Hastings Acceptance, using the joint density of position + momentum
    double top =  exp(joint_d(theta_prime, momentum_prime, log_pdf));
    double bottom =  exp(joint_d(chain.row(i-1), momentum, log_pdf));
    
    double alpha = top/bottom;
    
    if (R::runif(0,1) <= alpha){
      chain.row(i) = theta_prime;
      momentums.row(i) = momentum_prime;
      acceptances++;
    } else{
      chain.row(i) = chain.row(i-1);
      momentums.row(i) = momentums.row(i-1);
    }
  }
  
  return List::create(
    chain, 
    momentums, 
    ((double)(acceptances)/(double)(iterations))
  );
  
}

