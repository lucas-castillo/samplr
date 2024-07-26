#ifndef PDF_MANAGE
#define PDF_MANAGE

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

typedef std::function<double(Rcpp::NumericVector)> dfunc;

double abs_d(double x);

bool isClose(double a, double b, double tol=0.0001);

dfunc getPDF(
    const Rcpp::String &distr_name, 
    const Rcpp::List &distr_params, 
    const bool &log=false
);
  
double safe_log(const double &x);

dfunc getMixturePDF(
    std::vector<dfunc> &pdfs, 
    const Rcpp::NumericVector &weights, 
    const bool &logarithm = false
);

dfunc customPDF (const Rcpp::Function &f, const bool log = false);

dfunc managePDF(
    const Rcpp::StringVector &distr_name, 
    const Rcpp::List &distr_params, 
    const bool &isMix, 
    const Rcpp::NumericVector &weights, 
    const bool &log, 
    const Rcpp::Function &custom_func, 
    const bool &useCustom
);

Rcpp::NumericVector gridDensity_cpp(
    Rcpp::StringVector distr_name, Rcpp::List distr_params, bool isMix,
    Rcpp::NumericVector weights, Rcpp::NumericVector xxRange, 
    Rcpp::NumericVector yyRange, int cellsPerRow, Rcpp::Function densityFunc, 
    bool useCustomDensity
);

#endif