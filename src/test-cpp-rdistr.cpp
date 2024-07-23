// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include <testthat.h>
// 
#include "mh_utils.h"
#include "hmc_utils.h"
#include "pdf_manage.h"
#include "rdistr_manage.h"

#include "mh.h"
#include "mc3.h"
#include "hmc.h"
#include "mcrec.h"

using namespace Rcpp;

context("rdistr"){
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  
  test_that("get_rDistr"){
    // This distribution is not supported yet
    expect_error(get_rDistr("runif", List::create(1,3,2)));
  }
  
  
}
