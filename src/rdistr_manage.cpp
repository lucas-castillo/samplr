// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see custom_rDistr)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

// Random Generation Functions --------------

// Random Generation for build-in distributions
double get_rDistr(
    const String &distr_name, 
    const List &distr_params //, 
// const bool &log=false
)
{
  double sample_distr;
  // CONTINUOUS
  if (distr_name == "norm"){
    sample_distr = R::rnorm(distr_params(0), distr_params(1));
  } else {
    stop("This distribution is not supported yet.");
  }
  return sample_distr;
}


// Random Generation for custom distributions using rejection method
double custom_rDistr(
    const Function &f,
    const NumericVector &x_domains
){
  double sample_distr; // the sample from the distribution that will be returned
  NumericVector density(x_domains.length());
  bool reject = true;
  
  
  for (int i = 0; i < x_domains.length(); i++){
    try{
      density(i) = as<double>(f(x_domains(i)));
    } catch (const std::exception &exc){
      stop("\"x_domains\" contains values that are not defined in the distribution.");
    }
  }
  
  double max_dens = density(which_max(density));
  
  while(reject){
    int idx = floor(R::runif(0, x_domains.length()));
    double u = R::runif(0, 1);
    
    double density_value = as<double>(f(x_domains(idx)));
    if (!NumericVector::is_na(density_value) && density_value > 0 && u < density_value / max_dens) {
      sample_distr = x_domains(idx);
      reject = false;
    }
  }

  return sample_distr;
}


//[[Rcpp::export]]
double rDistr(
    const StringVector &distr_name,
    const List &distr_params,
    const Function &custom_func,
    const double &custom_start,
    const bool &useCustom
){
  double sample_distr;

  if (useCustom){
    sample_distr  = custom_start;
  } else {
    sample_distr = get_rDistr(distr_name(0), distr_params);
  }
  
  return sample_distr;
}


