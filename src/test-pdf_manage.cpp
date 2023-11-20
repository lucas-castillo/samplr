// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include <testthat.h>
//
#include "pdf_manage.h"


using namespace Rcpp;

context("PDF + Utils"){
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("isClose"){
    double x = 0;
    double y = .4;
    // Not close enough ...
    expect_false(isClose(x,y));
    // Except if we're really tolerant
    expect_true(isClose(x, y, 20));
  }
}
