#ifndef MH_UTILS
#define MH_UTILS

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// 
#include "pdf_manage.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

Rcpp::NumericVector alpha_trick(
    Rcpp::NumericVector random_jump,
    Rcpp::NumericVector last_jump,
    double alpha
);

Rcpp::NumericVector autocorrelated_metropolis_step_cpp(
    Rcpp::NumericMatrix &chain,
    Rcpp::NumericMatrix &proposals,
    Rcpp::NumericMatrix &jumps,
    Rcpp::NumericMatrix &true_jumps,
    const int &currentIndex,
    const double &last_prob,
    const Rcpp::NumericMatrix &sigma_prop,
    dfunc &pdf,
    const bool &discreteValues,
    const double &beta,
    const double &alpha
);

#endif