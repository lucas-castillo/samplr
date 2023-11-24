#ifndef HMC
#define HMC
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

Rcpp::List sampler_hmc_cpp(
    Rcpp::NumericVector start,
    Rcpp::StringVector distr_name,
    Rcpp::List distr_params,
    double epsilon,
    int L,
    int iterations,
    bool isMix,
    Rcpp::NumericVector weights,
    Rcpp::Function custom_func,
    bool useCustom
);

#endif
