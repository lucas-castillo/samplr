// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include "mc3.h"
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


using namespace Rcpp;


// ABS sampler//
//[[Rcpp::export]]
NumericVector subset_range(NumericVector x, int start, int end) {
  // Use the Range function to create a positional index sequence
  return x[Range(start, end)]; // end inclusive
}
//[[Rcpp::export]]
NumericVector support(NumericVector chain, double dec_bdry, int start_bias){
  NumericVector supportVector(chain.size()+1);
  supportVector(0) = start_bias;
  for (int i = 0; i < chain.size(); i++){
    if (chain(i) <= dec_bdry){
      supportVector(i+1) = -1;
    } else{
      supportVector(i+1) = 1;
    }
  }
  
  return(supportVector);
}
//[[Rcpp::export]]
NumericVector cumsum_sug(NumericVector x){
  NumericVector cumSupp = cumsum(x); // compute the cumulated number of evidence
  int n = cumSupp.size();
  return cumSupp[Range(1, n)];    // remove the first one as it is a start_bias
}
//[[Rcpp::export]]
NumericVector concatenate_vectors(NumericVector a, NumericVector b){
  if (a.size() == 0){
    return b;
  } else if (b.size() == 0){
    return a;
  }
  
  NumericVector c(a.size() + b.size());
  for (int i = 0; i<a.size(); i++){
    c(i) = a(i);
  }
  
  for (int i = 0; i<b.size(); i++){
    c(i + a.size()) = b(i);
  }
  return(c);
}
//[[Rcpp::export]]
int checkThreshold(NumericVector cumSupp, double delta, int caution){
  int supportPosition = -1;
  NumericVector cumSuppAbs = abs(cumSupp);
  for (int i = 0; i < cumSupp.size(); i++){
    if (cumSuppAbs(i) >= (delta + caution)){
      supportPosition = i;
      break;
    }
  }
  return supportPosition;
}

//[[Rcpp::export]]
NumericMatrix new_start_point(
    NumericVector chain, 
    int nChains, 
    int position, 
    int iterations
){
  NumericMatrix start_point(nChains, 1);
  
  for (int c = 0; c < nChains; c++){
    start_point(c, 0) = chain(position + c * iterations);
  }
  
  return(start_point);
}

// [[Rcpp::export]]
int give_me_sign(double x){
  if (x > 0){ 
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}

int bool_to_int(bool x){
  if (x){
    return 1;
  } else{
    return 0;
  }
}

// [[Rcpp::export]]
List ABS_sampler_cpp(
    NumericMatrix start_point,
    NumericVector trial_fdbk,
    StringVector distr_name, 
    int nChains,
    double dec_bdry,
    double d_sepn,
    double delta,
    double nd_time,
    double s_nd_time,
    double er_lambda,
    int mc3_iterations,
    double proposal_width
)
{ 
  int emp_ntrials = trial_fdbk.size();
  List sim_all_trials(emp_ntrials);
  int start_bias = 0;
  int caution = 0;
  Function f("rnorm"); // placeholder function
  
  int first_smpl_idx;
  List distr_params;
  List mc3_traces;
  NumericVector chain;
  NumericVector supportVector;
  NumericVector cumulative_support;
  NumericVector temp = {proposal_width};
  NumericMatrix sigma_prop(1, 1, temp.begin()); // 1x1 matrix with a proposal_width inside
  
  
  
  // Begin sampling
  for (int i = 0; i < emp_ntrials; i++){
    // Checking interruption every 1000 iterations
    if (i % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    NumericVector trialSamples;
    NumericVector trialSupport;
    int trialResponse;
    double trialNDTime;
    double trialDecisionTime;
    
    if (i == 0){
      first_smpl_idx = 0;
    } else {
      first_smpl_idx = 1;
    }
    
    if (trial_fdbk(i) == 1) {
      distr_params = List::create(d_sepn/2, 1);
    } else {
      distr_params = List::create(-1 * d_sepn/2, 1);
    }
    
    bool surpassedThreshold = false;
    while (!surpassedThreshold){
      mc3_traces = sampler_mc3_cpp(
        start_point, // start
        nChains, // nChains
        sigma_prop, // sigma_prop
        4, // delta_T
        true, // swap_all
        mc3_iterations, // iterations
        distr_name, // distr_name
        distr_params, // distr_params
        false, // discreteValues
        false, // isMix
        1, // weights
        f, // custom_func
        false // useCustom
      );
      
      chain = subset_range(mc3_traces[0], first_smpl_idx, mc3_iterations - 1); // cold chain
      
      supportVector = support(chain, dec_bdry, start_bias);
      cumulative_support = cumsum_sug(supportVector);
      int supportPosition = checkThreshold(cumulative_support, delta, caution);
      
      if (supportPosition == -1){
        start_point = new_start_point(
          mc3_traces[0],
                    nChains,
                    mc3_iterations - 1, // position = last Position
                    mc3_iterations
        );
        trialSamples = concatenate_vectors(trialSamples, chain);
        trialSupport = concatenate_vectors(trialSupport, cumulative_support);
        first_smpl_idx = 1;
        surpassedThreshold = false;
        // Rcout << "More iterations needed\n";
      } else{
        start_point = new_start_point(
          mc3_traces[0],
                    nChains,
                    supportPosition, // position = possition where boundary was crossed
                    mc3_iterations
        );
        
        trialSamples = concatenate_vectors(trialSamples, subset_range(chain, 0, supportPosition));
        trialSupport = concatenate_vectors(trialSupport, subset_range(cumulative_support, 0, supportPosition));
        
        surpassedThreshold = true;
        
      }
    }
    trialResponse = (give_me_sign(trialSupport(trialSupport.size() - 1)) + 1 )/ 2; // 0 for negative, 1 for positive
    trialNDTime = R::runif(nd_time, nd_time+s_nd_time);
    trialDecisionTime = R::rgamma(trialSamples.size(), 1/er_lambda);
    
    sim_all_trials(i) = List::create(
      _["trial"] = i + 1,
      _["samples"] = trialSamples,
      _["support"] = trialSupport,
      _["length"] = trialSamples.size(),
      _["response"] = trialResponse,
      _["feedback"] = trial_fdbk(i),
      _["accuracy"] = bool_to_int(trialResponse == trial_fdbk(i)),
      _["nd_time"] = trialNDTime,
      _["rt"] = trialDecisionTime + trialNDTime
    );
    start_bias = 0; // never start bias
    caution = 0; //1 - bool_to_int(trial_fdbk(i) == trialResponse);
  }
  return(sim_all_trials);
}

