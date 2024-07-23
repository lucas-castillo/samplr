// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see custom_rDistr)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

// Random Generation Functions --------------

// Random Generation for build-in distributions
double get_rDistr(
    const String &distr_name, 
    const List &distr_params //, 
// const bool &log=false
);

double rDistr(
    const StringVector &distr_name,
    const List &distr_params,
    const Function &custom_func,
    const double &custom_start,
    const bool &useCustom
);


