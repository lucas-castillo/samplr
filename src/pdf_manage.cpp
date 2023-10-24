// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;
typedef std::function<double(NumericVector)> dfunc;


dfunc getPDF(
    const String &distr_name, 
    const List &distr_params, 
    const bool &log=false
)
{
  dfunc pdf;
  // CONTINUOUS
  if (distr_name == "unif"){
    pdf = [distr_params, log](NumericVector x){return R::dunif(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "norm"){
    pdf = [distr_params, log](NumericVector x){return R::dnorm(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "lnorm"){
    pdf = [distr_params, log](NumericVector x){return R::dlnorm(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "gamma"){
    pdf = [distr_params, log](NumericVector x){return R::dgamma(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "beta"){
    pdf = [distr_params, log](NumericVector x){return R::dbeta(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "nbeta"){
    pdf = [distr_params, log](NumericVector x){return R::dnbeta(x(0), distr_params(0), distr_params(1), distr_params(2), log);};
  } else if (distr_name == "chisq"){
    pdf = [distr_params, log](NumericVector x){return R::dchisq(x(0), distr_params(0), log);};
  } else if (distr_name == "nchisq"){
    pdf = [distr_params, log](NumericVector x){return R::dnchisq(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "t"){
    pdf = [distr_params, log](NumericVector x){return R::dt(x(0), distr_params(0), log);};
  } else if (distr_name == "nt"){
    pdf = [distr_params, log](NumericVector x){return R::dnt(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "f"){
    pdf = [distr_params, log](NumericVector x){return R::df(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "nf"){
    pdf = [distr_params, log](NumericVector x){return R::dnf(x(0), distr_params(0), distr_params(1), distr_params(2), log);};
  } else if (distr_name == "cauchy"){
    pdf = [distr_params, log](NumericVector x){return R::dcauchy(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "exp"){
    pdf = [distr_params, log](NumericVector x){return R::dexp(x(0), distr_params(0), log);};
  } else if (distr_name == "logis"){
    pdf = [distr_params, log](NumericVector x){return R::dlogis(x(0), distr_params(0), distr_params(1), log);};
  } else if (distr_name == "weibull"){
    pdf = [distr_params, log](NumericVector x){return R::dweibull(x(0), distr_params(0), distr_params(1), log);};
    
    // RCPP-DIST Distributions
  } else if (distr_name == "4beta"){
    pdf = [distr_params, log](NumericVector x){return d4beta(x, distr_params(0), distr_params(1), distr_params(2), distr_params(3), log)[0];};
  } else if (distr_name == "lst"){
    pdf = [distr_params, log](NumericVector x){return dlst(x, distr_params(0), distr_params(1), distr_params(2), log)[0];};
  } else if (distr_name == "truncnorm"){
    pdf = [distr_params, log](NumericVector x){return dtruncnorm(x, distr_params(0), distr_params(1), distr_params(2), distr_params(3), log)[0];};
  } else if (distr_name == "trunct"){
    pdf = [distr_params, log](NumericVector x){return dtrunct(x, distr_params(0), distr_params(1), distr_params(2), log)[0];};
  } else if (distr_name == "trunclst"){
    pdf = [distr_params, log](NumericVector x){return dtrunclst(x, distr_params(0), distr_params(1), distr_params(2), distr_params(3), distr_params(4), log)[0];};
  } else if (distr_name == "triangular"){
    pdf = [distr_params, log](NumericVector x){return dtri(x, distr_params(0), distr_params(1), distr_params(2), log)[0];};
    // mv distributions use armadillo
  } else if (distr_name == "mvnorm"){
    pdf = [distr_params, log](NumericVector x){
      return dmvnorm(as<arma::rowvec>(x), as<arma::vec>(distr_params(0)), as<arma::mat>(distr_params(1)), log)[0];
    };
  } else if (distr_name == "mvt"){
    pdf = [distr_params, log](NumericVector x){
      return dmvt(as<arma::rowvec>(x), as<arma::vec>(distr_params(0)), as<arma::mat>(distr_params(1)), distr_params(2), log)[0];
    };
  }
  // DISCRETE
  else if (distr_name == "binom"){
    pdf = [distr_params, log](NumericVector x){return R::dbinom(x(0), distr_params(0), distr_params(1), log);};
  }else if (distr_name == "nbinom"){
    pdf = [distr_params, log](NumericVector x){return R::dnbinom(x(0), distr_params(0), distr_params(1), log);};
  }else if (distr_name == "nbinom_mu"){
    pdf = [distr_params, log](NumericVector x){return R::dnbinom_mu(x(0), distr_params(0), distr_params(1), log);};
  }else if (distr_name == "pois"){
    pdf = [distr_params, log](NumericVector x){return R::dpois(x(0), distr_params(0), log);};
  }else if (distr_name == "geom"){
    pdf = [distr_params, log](NumericVector x){return R::dgeom(x(0), distr_params(0), log);};
  }else if (distr_name == "hyper"){
    pdf = [distr_params, log](NumericVector x){return R::dhyper(x(0), distr_params(0), distr_params(1), distr_params(2), log);};
  }else if (distr_name == "wilcox"){
    pdf = [distr_params, log](NumericVector x){return R::dwilcox(x(0), distr_params(0), distr_params(1), log);};
  }else if (distr_name == "signrank"){
    pdf = [distr_params, log](NumericVector x){return R::dsignrank(x(0), distr_params(0), log);};
  }
  return pdf;
}

double safe_log(const double &x){
  if (x < 0){
    const unsigned int zero = 0;
    return log(zero);
  } else {
    return log(x);
  }
  
  
}

dfunc getMixturePDF(
    std::vector<dfunc> &pdfs, 
    const NumericVector &weights, 
    const bool &logarithm = false
){
  dfunc pdf;
  
  pdf = [pdfs, weights, logarithm](NumericVector x){
    double total_density = 0;
    for (unsigned i = 0; i < weights.size(); i++){
      dfunc p = pdfs[i];
      
      total_density +=  p(x) * weights[i];
    }
    if (logarithm){
      return safe_log(total_density);
    } else{
      return total_density;
    }
  };
  return pdf;
  
}

dfunc customPDF (const Function &f, const bool log = false){
  dfunc pdf = [f, log](NumericVector x){
    double d = 0;
    PutRNGstate();
    d = as<double>(f(x));
    GetRNGstate();
    if (log){
      d = safe_log(d);
    }
    return d;
  };
  return pdf;
}

dfunc managePDF(
    const StringVector &distr_name, 
    const List &distr_params, 
    const bool &isMix, 
    const NumericVector &weights, 
    const bool &log, 
    const Function &custom_func, 
    const bool &useCustom
){
  dfunc pdf;
  std::vector<dfunc> pdfs;
  
  if (useCustom){
    pdf  = customPDF(custom_func, log);
  } else if (!isMix){
    pdf = getPDF(distr_name(0), distr_params, log);
  } else {
    for (int i = 0; i < distr_name.size(); i++){
      pdfs.push_back(getPDF(distr_name(i), distr_params(i), false));
    }
    pdf = getMixturePDF(pdfs, weights, log);
  }
  
  return pdf;
}



// Plot Aid //
// [[Rcpp::export]]
NumericVector gridDensity_cpp(
    StringVector distr_name, List distr_params, bool isMix,
    NumericVector weights, NumericVector xxRange, NumericVector yyRange,
    int cellsPerRow, Function densityFunc, bool useCustomDensity
)
{
  dfunc pdf = managePDF(distr_name, distr_params, isMix, weights, false, densityFunc, useCustomDensity);
  
  NumericVector density(yyRange.size());
  for (int i = 0; i <yyRange.size(); i++){
    density(i) = pdf(NumericVector::create(xxRange[i],yyRange[i]));
  }
  return density;
}
