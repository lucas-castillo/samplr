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

bool isClose(double a, double b, double tol=0.0001){
  double diff = a - b;
  return abs(diff) < tol;
}
//[[Rcpp::export]]
NumericVector test_mh_step(){
    NumericMatrix chain(2, 1);
    NumericMatrix proposals(2, 1);
    NumericMatrix jumps(2, 1);
    NumericMatrix true_jumps(2, 1);
    int currentIndex = 1;
    double last_prob = .5;
    NumericMatrix sigma_prop(1);
    dfunc pdf = getPDF("norm", List::create(0, 1), false);
    bool discreteValues = false;
    double beta = 0;
    double alpha = 0;

    return autocorrelated_metropolis_step_cpp(
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

}



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
    NumericVector zero = {0};
    LogicalVector same = res==zero;
    print(res);
    expect_true(same(0));
    
  }
  // test_that("autocorrelated metropolis step"){
  //   NumericMatrix chain(2, 1);
  //   NumericMatrix proposals(2, 1);
  //   NumericMatrix jumps(2, 1);
  //   NumericMatrix true_jumps(2, 1);
  //   int currentIndex = 1;
  //   double last_prob = .5;
  //   NumericMatrix sigma_prop(1);
  //   dfunc pdf = getPDF("norm", List::create(0, 1), false);
  //   bool discreteValues = false;
  //   double beta = 0;
  //   double alpha = 0;
  //   
  //   autocorrelated_metropolis_step_cpp(
  //       chain,
  //       proposals,
  //       jumps,
  //       true_jumps,
  //       currentIndex,
  //       last_prob,
  //       sigma_prop,
  //       pdf,
  //       discreteValues,
  //       beta,
  //       alpha
  //   )
  //   
  // }
  
}