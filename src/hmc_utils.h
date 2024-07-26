#ifndef HMC_UTILS
#define HMC_UTILS

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// 
#include "pdf_manage.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


Rcpp::NumericVector gradient(dfunc &func, const Rcpp::NumericVector &x, double Temp = 1);
double dotProduct(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y);
double joint_d(
    const Rcpp::NumericVector &theta, 
    const Rcpp::NumericVector &momentum, 
    dfunc &log_func, 
    double Temp=1
);
void leapfrog_step_cpp(
    Rcpp::NumericVector &theta, 
    Rcpp::NumericVector &momentum, 
    const double &epsilon, 
    dfunc &log_pdf, 
    const int &L, 
    double Temp=1
);
Rcpp::NumericVector RecycledMomentumUpdate(
    Rcpp::NumericVector momentum,
    double alpha,
    double Temp=1
);
Rcpp::NumericVector drawMomentum(int dimensions);
#endif