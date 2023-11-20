// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include <testthat.h>
// 
#include "mh_utils.h"
#include "pdf_manage.h"


using namespace Rcpp;


context("MH Utils"){
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("alpha update") {
    double a = .5;
    NumericVector l = {2};
    NumericVector r = -l / pow(3, .5);
    NumericVector res = alpha_trick(r,l,a);
    expect_true(isClose(res(0), 0));
    
  }
  test_that("autocorrelated metropolis step"){
    NumericVector v = {1, 0};
    NumericVector v2 = {1};
    NumericMatrix chain(2, 1, v.begin());
    NumericMatrix proposals(2, 1, v.begin());
    NumericMatrix jumps(2, 1, v.begin());
    NumericMatrix true_jumps(2, 1, v.begin());
    int currentIndex = 1;
    NumericMatrix sigma_prop(1, 1, v2.begin());
    dfunc pdf = getPDF("norm", List::create(0, 1), false);
    double last_prob = .5;
    bool discreteValues = false;
    double beta = 0;
    double alpha = 1;
    
    NumericVector result = autocorrelated_metropolis_step_cpp(
      chain,
      proposals,
      jumps,
      true_jumps,
      currentIndex,
      last_prob,
      sigma_prop,
      pdf,
      discreteValues,
      beta,
      alpha
    );
    // With Alpha = 1, the jump will be identical to last (1) -- so returned pdf is dnorm(2)
    expect_true(isClose(result(0),R::dnorm(2, 0, 1, false)));
    
    // With beta = 0, the jump will be always accepted
    expect_true(isClose(result(1), 1));
    
    beta = 10;
    last_prob = 1000;
    
    result = autocorrelated_metropolis_step_cpp(
      chain,
      proposals,
      jumps,
      true_jumps,
      currentIndex,
      last_prob,
      sigma_prop,
      pdf,
      discreteValues,
      beta,
      alpha
    );
    
    // With beta = 1 and last_prob = 10 the result is always rejected.
    expect_true(isClose(result(1), 0));
    
    
    // And this returns the last_prob we provided
    expect_true(isClose(result(0), last_prob));
    
    last_prob = -10;
    result = autocorrelated_metropolis_step_cpp(
      chain,
      proposals,
      jumps,
      true_jumps,
      currentIndex,
      last_prob,
      sigma_prop,
      pdf,
      discreteValues,
      beta,
      alpha
    );
    
    // A negative "last prob" will return the proposal, always
    expect_true(isClose(result(1), 1));
    
    NumericVector v3 = {1.2, 0};
    NumericMatrix jumps2(2, 1, v3.begin()); // jump = 1.2
    discreteValues = true; // this will round the proposal
    
    result = autocorrelated_metropolis_step_cpp(
      chain,
      proposals,
      jumps,
      true_jumps,
      currentIndex,
      last_prob,
      sigma_prop,
      pdf,
      discreteValues,
      beta,
      alpha
    );
    // The pdf will be after jump of +1, not +1.2    
    expect_true(isClose(result(0),R::dnorm(2, 0, 1, false)));
    
  }
  
}
