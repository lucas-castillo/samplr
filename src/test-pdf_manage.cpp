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
  
  test_that("getPDF"){
    dfunc pdf;
    // Univariate Continuous
    pdf = getPDF("unif", List::create(0,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 1));
    pdf = getPDF("norm", List::create(0,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.3989423));
    pdf = getPDF("lnorm", List::create(0,1), false);
    expect_true(isClose(pdf(NumericVector::create(1)), 0.3989423));
    pdf = getPDF("gamma", List::create(1,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 1));
    pdf = getPDF("beta", List::create(2,2), false);
    expect_true(isClose(pdf(NumericVector::create(.5)), 1.5));
    pdf = getPDF("nbeta", List::create(2,2,2), false);
    expect_true(isClose(pdf(NumericVector::create(.5)), 1.402602));
    pdf = getPDF("chisq", List::create(2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), .5));
    pdf = getPDF("nchisq", List::create(2, 2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.1839397));
    pdf = getPDF("t", List::create(2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.3535534));
    pdf = getPDF("nt", List::create(2,2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.04784825));
    pdf = getPDF("f", List::create(2,4), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 1));
    pdf = getPDF("nf", List::create(2,4,2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.3678794));
    pdf = getPDF("cauchy", List::create(0,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.3183099));
    pdf = getPDF("exp", List::create(1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 1));
    pdf = getPDF("logis", List::create(0,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), .25));
    pdf = getPDF("weibull", List::create(1,1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 1));
    pdf = getPDF("4beta", List::create(3,3,2,6), false);
    expect_true(isClose(pdf(NumericVector::create(3)), 0.2636719));
    pdf = getPDF("lst", List::create(2,0, 1), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.3535534));
    pdf = getPDF("truncnorm", List::create(0,1,-2, 2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.4179596));
    pdf = getPDF("trunct", List::create(2, -2, 2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.4330127));
    pdf = getPDF("trunclst", List::create(2, 0,1, -2,2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.4330127));
    pdf = getPDF("triangular", List::create(1,3,2), false);
    expect_true(isClose(pdf(NumericVector::create(1.5)), .5));
    
    // Multivariate Continuous
    NumericVector zeroes = {0,0};
    NumericVector cov_v = {1, 0, 0, 1};
    NumericMatrix cov( 2 , 2 , cov_v.begin() );
    
    pdf = getPDF("mvnorm", List::create(zeroes, cov), false);
    expect_true(isClose(pdf(zeroes), 0.1591549));
    
    pdf = getPDF("mvt", List::create(zeroes, cov, 3), false);
    expect_true(isClose(pdf(zeroes), 0.1591549));
    
    // Univariate Discrete
    pdf = getPDF("binom", List::create(3,.5), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.125));
    pdf = getPDF("nbinom", List::create(3,.5), false);
    expect_true(isClose(pdf(NumericVector::create(5)), 0.08203125));
    pdf = getPDF("nbinom_mu", List::create(4,4), false);
    expect_true(isClose(pdf(NumericVector::create(5)), 0.109375));
    pdf = getPDF("pois", List::create(2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.1353353));
    pdf = getPDF("geom", List::create(.5), false);
    expect_true(isClose(pdf(NumericVector::create(2)), .125));
    pdf = getPDF("hyper", List::create(3,4,5), false);
    expect_true(isClose(pdf(NumericVector::create(2)), 0.5714286));
    pdf = getPDF("wilcox", List::create(2,3), false);
    expect_true(isClose(pdf(NumericVector::create(2)), .2));
    pdf = getPDF("signrank", List::create(2), false);
    expect_true(isClose(pdf(NumericVector::create(0)), .25));
  }
  
  test_that("safe_log"){
    double neg_inf = safe_log(-2);
    expect_true(neg_inf == R_NegInf);
    expect_true(safe_log(2) == log(2));
  }
  
  test_that("getMixturePDF"){
    dfunc pdf1 = getPDF("norm", List::create(-1,1), false);
    dfunc pdf2 = getPDF("norm", List::create(0,2), false);
    std::vector<dfunc> v{pdf1, pdf2};
    
    dfunc pdf = getMixturePDF(v, NumericVector::create(.6, .4), false);
    expect_true(isClose(pdf(NumericVector::create(0)), 0.2249709));
    
    dfunc pdf_log = getMixturePDF(v, NumericVector::create(.6, .4), true);
    expect_true(isClose(pdf_log(NumericVector::create(0)), -1.491784));
    
  }
}
