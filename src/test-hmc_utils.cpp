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
    dfunc log_pdf = getPDF("norm", List::create(0, 1), true);
    expect_true(isClose(gradient(log_pdf, NumericVector::create(0))(0), 0));
    
    // Expected gradient is -x
    expect_true(isClose(gradient(log_pdf, NumericVector::create(1))(0), -1));
    
    // Gradient decreases if temperature increases
    expect_true(isClose(gradient(log_pdf, NumericVector::create(1), 2)(0), -.5));
  }
  
  test_that("dotProduct"){
    NumericVector x = {1,2};
    NumericVector y = {3,4};
    NumericVector z = {5};
    
    expect_true(dotProduct(x,y) == 11);
    expect_error(dotProduct(x,z));
  }
  
  test_that("joint_d"){
    dfunc log_pdf = getPDF("norm", List::create(0, 1), true);
    NumericVector x = {0};
    NumericVector p = {3};
    
    expect_true(isClose(joint_d(x, p, log_pdf), -5.418939));
    expect_true(isClose(joint_d(x, p, log_pdf, 3), -0.8063128));
    
  }
  
  test_that("leapfrog_step_cpp"){
    NumericVector theta = {0};
    NumericVector momentum = {1};
    double epsilon = .1;
    dfunc log_pdf = getPDF("norm", List::create(0, 1), true);
    int L = 10;
    leapfrog_step_cpp(theta, momentum, epsilon, log_pdf, L);
    expect_true(isClose(theta(0), 0.8427504));
    expect_true(isClose(momentum(0), -0.5399513));
    // this is reversible
    leapfrog_step_cpp(theta, momentum, epsilon, log_pdf, L);
    expect_true(isClose(theta(0), 0));
    expect_true(isClose(momentum(0), 1));

    // Temperature
    leapfrog_step_cpp(theta, momentum, epsilon, log_pdf, L, 2);
    expect_true(isClose(theta(0), 0.9194587));
    expect_true(isClose(momentum(0), -0.7601488));
  }
 }
