// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include <testthat.h>
// 
#include "hmc_utils.h"
#include "pdf_manage.h"


using namespace Rcpp;

context("HMC Utils"){
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Gradient") {
    // Expected gradient = 0
    dfunc pdf = getPDF("norm", List::create(0, 1), false);
    expect_true(isClose(gradient(pdf, NumericVector::create(0))(0), 0));
    
    // Expected gradient approx -0.24
    double pi = atan(1)*4;
    double g = - 1 / (sqrt(2 * exp(1) * pi));
    expect_true(isClose(gradient(pdf, NumericVector::create(1))(0), g));
    
    // Gradient decreases if temperature increases
    expect_true(isClose(gradient(pdf, NumericVector::create(1), 2)(0), g/2));
  }
 }
