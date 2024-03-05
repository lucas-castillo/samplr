// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include "mc3.h"
#include "mcrec.h"
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


using namespace Rcpp;


// ABS sampler//
NumericVector subset_range(NumericVector x, int start, int end) {
  // Use the Range function to create a positional index sequence
  return x[Range(start, end)]; // end inclusive
}

NumericVector support(NumericVector chain, double dec_bdry, int prior_bias){
  NumericVector supportVector(chain.size()+1);
  supportVector(0) = prior_bias;
  for (int i = 0; i < chain.size(); i++){
    if (chain(i) <= dec_bdry){
      supportVector(i+1) = -1;
    } else{
      supportVector(i+1) = 1;
    }
  }
  
  return(supportVector);
}

NumericVector cumsum_sug(NumericVector x){
  NumericVector cumSupp = cumsum(x); // compute the cumulated number of evidence
  int n = cumSupp.size();
  return cumSupp[Range(1, n)];    // remove the first one as it is a prior_bias
}

NumericVector concatenate_vectors(NumericVector a, NumericVector b){
  if (a.size() == 0){
    return b;
  } else if (b.size() == 0){
    return a;
  }
  
  NumericVector c(a.size() + b.size());
  for (int i = 0; i<a.size(); i++){
    c(i) = a(i);
  }
  
  for (int i = 0; i<b.size(); i++){
    c(i + a.size()) = b(i);
  }
  return(c);
}

int checkThreshold(NumericVector cumSupp, double delta){
  int supportPosition = -1;
  NumericVector cumSuppAbs = abs(cumSupp);
  for (int i = 0; i < cumSupp.size(); i++){
    if (cumSuppAbs(i) >= (delta)){
      supportPosition = i;
      break;
    }
  }
  return supportPosition;
}


NumericMatrix new_start_point(
    NumericVector chain, 
    int n_chains, 
    int position, 
    int iterations
){
  NumericMatrix start_point(n_chains, 1);
  
  for (int c = 0; c < n_chains; c++){
    start_point(c, 0) = chain(position + c * iterations);
  }
  
  return(start_point);
}


int give_me_sign(double x){
  if (x > 0){ 
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}

int bool_to_int(bool x){
  if (x){
    return 1;
  } else{
    return 0;
  }
}

