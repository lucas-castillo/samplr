// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>

// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>
#include "mc3.h"
#include "mcrec.h"
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


using namespace Rcpp;


// ABS sampler//
NumericVector subset_range(NumericVector x, int start, int end) {
  // Use the Range function to create a positional index sequence
  return x[Range(start, end)]; // end inclusive
}

NumericVector support(NumericVector chain, double dec_bdry, int prior_bias){
  NumericVector supportVector(chain.size()+1);
  supportVector(0) = prior_bias;
  for (int i = 0; i < chain.size(); i++){
    if (chain(i) <= dec_bdry){
      supportVector(i+1) = -1;
    } else{
      supportVector(i+1) = 1;
    }
  }
  
  return(supportVector);
}

NumericVector cumsum_sug(NumericVector x){
  NumericVector cumSupp = cumsum(x); // compute the cumulated number of evidence
  int n = cumSupp.size();
  return cumSupp[Range(1, n)];    // remove the first one as it is a prior_bias
}

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

int checkThreshold(NumericVector cumSupp, double delta){
  int supportPosition = -1;
  NumericVector cumSuppAbs = abs(cumSupp);
  for (int i = 0; i < cumSupp.size(); i++){
    if (cumSuppAbs(i) >= (delta)){
      supportPosition = i;
      break;
    }
  }
  return supportPosition;
}


NumericMatrix new_start_point(
    NumericVector chain, 
    int n_chains, 
    int position, 
    int iterations
){
  NumericMatrix start_point(n_chains, 1);
  
  for (int c = 0; c < n_chains; c++){
    start_point(c, 0) = chain(position + c * iterations);
  }
  
  return(start_point);
}


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

//[[Rcpp::export]]
List Zhu23ABS_cpp(
    int task_id, // 1 represent estimate task, 2 represent two-alternative choice task
    NumericVector trial_stim,
    NumericMatrix start_point,
    StringVector distr_name,
    double proposal_width,
    int n_chains,
    NumericVector prior_on_resp,
    int stop_rule, 
    double nd_time, 
    double s_nd_time,
    double er_lambda,
    int mc3_iterations = 100, // these three parameters are unused in task 1
    double dec_bdry = 0,
    double discrim = 0
){
  
  int emp_ntrials = trial_stim.size();
  List sim_all_trials(emp_ntrials);
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
    int prior_bias;
    double trialNDTime;
    double trialDecisionTime;
    
    // set the index of the first sample for a MC3 chain
    if (i == 0){
      first_smpl_idx = 0;
    } else {
      first_smpl_idx = 1;
    }
    
    switch (task_id){
    case 1:
      distr_params = List::create(trial_stim(i), 4);
      
      mc3_traces = sampler_mc3_cpp(
        start_point, // start
        n_chains, // n_chains
        sigma_prop, // sigma_prop
        4, // delta_T
        true, // swap_all
        stop_rule + first_smpl_idx, // iterations
        distr_name, // distr_name
        distr_params, // distr_params
        false, // discreteValues
        false, // isMix
        1, // weights
        f, // custom_func
        false // useCustom
      );
      
      start_point = new_start_point(
        mc3_traces[0],
                  n_chains,
                  stop_rule + first_smpl_idx - 1, // position = position where boundary was crossed
                  stop_rule + first_smpl_idx
      );
      
      trialSamples = subset_range(mc3_traces[0], first_smpl_idx, stop_rule + first_smpl_idx - 1); // cold chain
      trialNDTime = R::runif(nd_time, nd_time+s_nd_time);
      trialDecisionTime = R::rgamma(stop_rule, 1/er_lambda);
      
      sim_all_trials(i) = List::create(
        _["trial"] = i + 1,
        _["samples"] = trialSamples,
        _["stimulus"] = trial_stim(i),
        _["nd_time"] = trialNDTime,
        _["rt"] = trialDecisionTime + trialNDTime
      );
      
      break;
      
    case 2: // force-choice task
      // for the force-choice task, the trial_stim(i) should be either 1 or 2
      if (trial_stim(i) == 1) {
        distr_params = List::create(-1 * discrim/2, 1);
      } else if (trial_stim(i) == 2) {
        distr_params = List::create(discrim/2, 1);
      } else {
        stop("trial_stim has more than two levels.");
      }
      
      bool surpassedThreshold = false;
      while (!surpassedThreshold){
        mc3_traces = sampler_mc3_cpp(
          start_point, // start
          n_chains, // n_chains
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
        
        //Rcout << "The prior is" << prior_on_resp << "\n";
        prior_bias = prior_on_resp(1) - prior_on_resp(0);
        
        supportVector = support(chain, dec_bdry, prior_bias);
        cumulative_support = cumsum_sug(supportVector);
        int supportPosition = checkThreshold(cumulative_support, stop_rule);
        
        if (supportPosition == -1){
          start_point = new_start_point(
            mc3_traces[0],
                      n_chains,
                      mc3_iterations - 1, // position = last Position
                      mc3_iterations
          );
          trialSamples = concatenate_vectors(trialSamples, chain);
          trialSupport = concatenate_vectors(trialSupport, cumulative_support);
          first_smpl_idx = 1;
          
          if (trialSamples.size() < 999){
            surpassedThreshold = false;
          }else{
            stop("The simulated sequence is two long. Please examine your parameters.");
          }
          
          // Rcout << "More iterations needed\n";
        } else{
          start_point = new_start_point(
            mc3_traces[0],
                      n_chains,
                      supportPosition, // position = position where boundary was crossed
                      mc3_iterations
          );
          
          trialSamples = concatenate_vectors(trialSamples, subset_range(chain, 0, supportPosition));
          trialSupport = concatenate_vectors(trialSupport, subset_range(cumulative_support, 0, supportPosition));
          
          surpassedThreshold = true;
        }
      }
      
      double trialLength = trialSamples.size();
      trialResponse = (give_me_sign(trialSupport(trialLength - 1)) + 3 )/ 2; // 1 for lower than the dec_bdry, 2 for higher than the dec_bdry
      trialNDTime = R::runif(nd_time, nd_time+s_nd_time);
      trialDecisionTime = R::rgamma(trialLength, 1/er_lambda);
      
      double evidDifference = trialSupport(trialLength - 1);
      double conf_posi = (prior_on_resp[1] + (trialLength + evidDifference - prior_bias) / 2)/(trialLength + prior_on_resp[0] + prior_on_resp[1]);// The confidence of "positive" response
      
      // update the prior on responses
      prior_on_resp = NumericVector::create(1, 1); // reset
      prior_on_resp(trial_stim(i)-1) = prior_on_resp(trial_stim(i)-1) + 1; // update the prior on responses based on the stimulus
      
      sim_all_trials(i) = List::create(
        _["trial"] = i + 1,
        _["samples"] = trialSamples,
        //_["support"] = trialSupport,
        //_["length"] = trialLength,
        _["response"] = trialResponse,
        _["stimulus"] = trial_stim(i),
        _["accuracy"] = bool_to_int(trialResponse == trial_stim(i)),
        //_["nd_time"] = trialNDTime,
        _["rt"] = trialDecisionTime + trialNDTime,
        _["confidence"] = std::max(conf_posi, 1-conf_posi)
      );
      break;
    }
    
  }
  
  return(sim_all_trials);
}


