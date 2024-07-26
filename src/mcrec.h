#ifndef MCREC
#define MCREC
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


#include "hmc_utils.h"
#include "pdf_manage.h"

Rcpp::List sampler_mc_rec_cpp(
    Rcpp::NumericVector start,
    int nChains, // for non-tempered versions of HMC and HOR, set this to one
    double delta_T, // for all chains having the same temp, set to 0
    bool swap_all,
    double iterations,
    Rcpp::StringVector distr_name,
    Rcpp::List distr_params,
    bool discreteValues,
    bool isMix,
    Rcpp::NumericVector weights,
    Rcpp::Function custom_func,
    bool useCustom,
    double epsilon,
    int L,
    double alpha // for HMC, set this to zero
);

#endif
