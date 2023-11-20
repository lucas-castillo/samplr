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

#include "mh.h"
#include "mc3.h"
#include "hmc.h"
#include "mcrec.h"

using namespace Rcpp;

context("CPP Samplers"){
  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("MH"){
    Function f("rnorm");
    NumericVector one = {1};
    NumericMatrix one_matrix(1, 1, one.begin());
    List res = sampler_mh_cpp(
      NumericVector::create(0), //Rcpp::NumericVector start,
      one_matrix, //Rcpp::NumericMatrix sigma_prop,
      100, //int iterations,
      StringVector::create("norm"), //Rcpp::StringVector distr_name,
      List::create(0,1), //Rcpp::List distr_params,
      false, //bool discreteValues,
      false, //bool isMix,
      NumericVector::create(1), //Rcpp::NumericVector weights,
      f, //Rcpp::Function custom_func,
      false, //bool useCustom,
      .5 //double alpha=0
    );
  }
  
}
