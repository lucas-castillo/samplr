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

int bool_to_int(bool x){
  if (x){return 1;}
  else {return 0;}
}

// convert a double variable x to a n_rep * 1 matrix
NumericMatrix double_to_matrix(double x, int n_rep){
  NumericMatrix m(n_rep, 1);
  for (int i = 0; i < n_rep; i++){
    m(i, 0) = x;
  }
  return(m);
}

// A function that returns one random value with a given distribution and related parameters
double rDistr(StringVector distr_name, List distr_params) {
  if (distr_name(0) == "norm") {
    return R::rnorm(distr_params(0), distr_params(1));
  } else {
    stop("Distribution not supported.");
  }
}

// A function for getting the last sample from the chain
NumericMatrix mc3_last_sample(NumericMatrix chain, int stop_position, int nChains, int total_iterations){
  NumericMatrix last_sample(nChains, 1);
  
  for (int i=0; i<nChains; i++){
    last_sample(i, 0) = chain(stop_position + i*total_iterations, 0);
  }
  return(last_sample);
}



// Samplers for the Zhu23ABS with the relative stopping rule
List sampler_mc3_rltv_stop_cpp(
    NumericMatrix start, // Numeric Matrix of starts for each chain.
    int nChains,
    NumericMatrix sigma_prop,
    double delta_T,
    bool swap_all,
    double iterations,
    StringVector distr_name,
    List distr_params,
    bool discreteValues,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom,
    int stop_rule, // these four parameters are for the ABS
    int first_sample_idx,
    NumericVector acc_evid,
    double dec_bdry
){
  // init the vars of mc3 ------------------------------------------------
  NumericVector acceptances(nChains);
  double alpha=0;
  int swap_attempts = 0;
  int swap_accepts= 0;
  NumericMatrix swaps(iterations*nChains, 3);
  int n_dim = start.ncol();
  
  NumericMatrix chain(iterations*nChains, n_dim);
  NumericMatrix proposals(iterations*nChains, n_dim);
  NumericVector beta(nChains);
  NumericMatrix ps(nChains, iterations);
  NumericMatrix jumpsHistory(iterations*nChains, n_dim);
  NumericMatrix jumpsUse(iterations*nChains, n_dim);
  NumericVector momentum;
  
  dfunc pdf = managePDF(distr_name, distr_params, isMix, weights, false, custom_func, useCustom);
  for (int i = 0; i < nChains; i++){
    chain.row(0 + iterations * i) = start.row(i);
    beta(i) = 1 / (1 + delta_T * i);
    ps(i,0) = pdf(start.row(i));
  }
  
  int nSwaps;
  if (swap_all){
    nSwaps = floor(nChains/2);
  } else {
    nSwaps = 1;
  }
  
  NumericVector v(nChains);
  v = seq(0, nChains - 1);
  
  // init the vars of ABS ------------------------------------------------
  int trialResponse;
  NumericVector trialSamples(iterations); // the cold chain
  bool start_sampling = true;
  
  trialSamples(0) = start(0, 0);
  
  // evaluate the start point if it is considered as one of the samples
  if (first_sample_idx == 0){
    if (trialSamples(0) <= dec_bdry) {
      acc_evid(0) += 1;
    } else {
      acc_evid(1) += 1;
    }
  }
  
  if ( acc_evid(0) - acc_evid(1) == stop_rule ){
    trialResponse = 1;
    trialSamples = trialSamples(0);
    start_sampling = false;
  } else if (acc_evid(1) - acc_evid(0) == stop_rule) {
    trialResponse = 2;
    trialSamples = trialSamples(0);
    start_sampling = false;
  }
  
  // Start sampling
  if (start_sampling){
    // run the sampler ------------------------------------------------
    for (int i = 1; i < iterations; i++){
      
      for (int ch = 0; ch < nChains; ch++){
        NumericVector accept;
        if (i == 1){
          accept =  autocorrelated_metropolis_step_cpp(
            chain, // NumericMatrix &chain, 
            proposals, // NumericMatrix &proposals, 
            jumpsUse, // NumericMatrix &jumps, 
            jumpsHistory,// NumericMatrix &true_jumps, 
            i + iterations * ch, // const int &currentIndex, 
            ps(ch, i-1), // const double &lastP, 
            sigma_prop, // const NumericMatrix &sigma_prop, 
            pdf, // dfunc &pdf, 
            discreteValues, // const bool &discreteValues, 
            beta(ch), // const double &beta, 
            0  // const double &alpha
          );  
        } else{
          accept =  autocorrelated_metropolis_step_cpp(
            chain, // NumericMatrix &chain, 
            proposals, // NumericMatrix &proposals, 
            jumpsUse, // NumericMatrix &jumps, 
            jumpsHistory,// NumericMatrix &true_jumps, 
            i + iterations * ch, // const int &currentIndex, 
            ps(ch, i-1), // const double &lastP, 
            sigma_prop, // const NumericMatrix &sigma_prop, 
            pdf, // dfunc &pdf, 
            discreteValues, // const bool &discreteValues, 
            beta(ch), // const double &beta, 
            alpha // const double &alpha
          );  
        }
        
        ps(ch, i) = accept(0);
        acceptances(ch) = acceptances(ch) + accept(1);
      }
      
      // once a step in every chain has been done, proceed to swap chains
      if (nChains > 1){
        // arrange chains randomly
        
        v = sample(v, nChains, false);
        
        // swap nSwaps times (depending on swap_all)
        for (int k = 0; k < nSwaps; k++){
          swap_attempts++;
          int m = v[k*2];
          int n = v[k*2 + 1];
          // chains are swapped with probability alpha, which is the ratio between:
          // - the product of the density of each chain's location at the temperature of the other chain, and
          // - the product of the density of each chain's location at their own temperatures
          NumericVector t1 = chain.row(i+iterations * m);
          NumericVector t2 = chain.row(i+iterations * n);
          double m_pdf = pdf(chain.row(i + iterations * m));
          double n_pdf = pdf(chain.row(i + iterations * n));
          
          double top = pow(m_pdf, beta(n)) * pow(n_pdf, beta(m));
          double bottom = pow(m_pdf,beta(m)) * pow(n_pdf, beta(n));
          
          
          if ((bottom != 0 && R::runif(0,1) <= top/bottom) || (bottom == 0 && top > 0)){
            
            // Swap Positions
            NumericVector temp = chain.row(i + iterations * m);
            chain.row(i + iterations * m) = chain.row(i + iterations * n);
            chain.row(i + iterations * n) = temp;
            // Record Swap
            swaps.row(swap_accepts) = NumericVector::create(i+1, m, n);
            swap_accepts++;
            // Swap probabilities
            t1 = chain.row(i+iterations * m);
            t2 = chain.row(i+iterations * n);
            double TEMP = ps(m, i);
            ps(m,i) = ps(n, i);
            ps(n,i) = TEMP;
            // Swap last jump
            if (i != (iterations-1)){
              // Reset jump after swap
              NumericVector zeros (n_dim);
              arma::mat random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
              NumericVector random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
              
              jumpsUse.row(i + iterations * m) = random_jump;
              
              
              random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
              random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
              
              jumpsUse.row(i + iterations * n) = random_jump;
            }
          }
        }
      }
      
      // evaluate the new sample
      double new_sample = chain(i, 0);
      trialSamples(i) = new_sample;
      
      if (new_sample <= dec_bdry) {
        acc_evid(0) += 1;
      } else {
        acc_evid(1) += 1;
      }
      
      if ( acc_evid(0) - acc_evid(1) == stop_rule ){
        trialResponse = 1;
        trialSamples = trialSamples[Range(first_sample_idx, i)];
        break; // break the loop
      } else if (acc_evid(1) - acc_evid(0) == stop_rule) {
        trialResponse = 2;
        trialSamples = trialSamples[Range(first_sample_idx, i)];
        break;
      }
      
      // if the sampling process reaches to the maximum limitation
      if (i == iterations-1){
        warning("The simulated sequence is two long. Please examine your parameters.");
        trialResponse = R_NaN;
      }
      
    } // the end of sampling
    
  }// the end of if (start_sampling)
  
  return List::create(
    _["chain"] = chain,
    _["samples"] = trialSamples,
    _["response"] = trialResponse,
    _["acc_evid"] = acc_evid
  );
}



//[[Rcpp::export]]
List Zhu23ABS_cpp(
    int task_id, // 1 represent the fixed stopping rule, 2 represent the relative stopping rule
    NumericVector trial_stim,
    StringVector distr_name,
    NumericVector distr_add_params,
    double proposal_width,
    int n_chains,
    NumericVector provided_start_point,
    int stop_rule, 
    double nd_time, 
    double s_nd_time,
    double lambda,
    NumericVector prior_on_resp = NumericVector::create(1,1), // the following five parameters are only for the relative stopping rule
    bool prior_depend = true,
    int mc3_iterations = 1000,
    double dec_bdry = 0,
    double discrim = 0
){
  
  int emp_ntrials = trial_stim.size();
  List sim_all_trials(emp_ntrials);
  Function f("rnorm"); // placeholder function
  NumericVector trialSamples(1);
  List sampler_results;
  NumericMatrix start_point_m;
  NumericVector acc_evid = clone(prior_on_resp); // accumulated evidence
  NumericMatrix sigma_prop = double_to_matrix(proposal_width, 1);// 1x1 matrix with a proposal_width inside
  
  // Begin sampling
  for (int i = 0; i < emp_ntrials; i++){
    // Checking interruption every 1000 iterations
    if (i % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    List distr_params;
    double start_point;
    int trialResponse;
    double trialNDTime;
    double trialDecisionTime;
    double conf;
    int first_sample_idx; // either 0 or 1, 0 means the start point is the first sample,
    
    switch (task_id){
    case 1:  // fixed stopping rule -----------------------------------------------------------------------------------
      
      distr_params = List::create(trial_stim(i), distr_add_params(i));
      
      // set the start point
      
      if (all(is_na(provided_start_point))){ // if users did not provide any start points
        if (i == 0){
          start_point = rDistr(distr_name, distr_params);
          start_point_m = double_to_matrix(start_point, n_chains); // convert start_point to a matrix
          first_sample_idx = 0;
        } else {
          first_sample_idx = 1;
        }
      } else {
        start_point = provided_start_point(i);
        start_point_m = double_to_matrix(start_point, n_chains);
        first_sample_idx = 0;
      }
      
      sampler_results = sampler_mc3_cpp(
        start_point_m, // start
        n_chains, // n_chains
        sigma_prop, // sigma_prop
        4, // delta_T
        true, // swap_all
        stop_rule + first_sample_idx, // iterations
        distr_name, // distr_name
        distr_params, // distr_params
        false, // discreteValues
        false, // isMix
        1, // weights
        f, // custom_func
        false // useCustom
      );
      
      trialSamples = subset_range(sampler_results["chain"], first_sample_idx, stop_rule + first_sample_idx - 1); // cold chain
      start_point_m = mc3_last_sample(sampler_results["chain"], stop_rule + first_sample_idx - 1, n_chains, stop_rule + first_sample_idx); // the start_point for the next trial
      trialNDTime = R::runif(nd_time, nd_time+s_nd_time);
      trialDecisionTime = R::rgamma(stop_rule, 1/lambda);
      
      sim_all_trials(i) = List::create(
        _["trial"] = i + 1,
        _["samples"] = trialSamples,
        _["stimulus"] = trial_stim(i),
        _["rt"] = trialDecisionTime + trialNDTime
      );
      
      
      break;
      
    case 2: // relative stopping rule ----------------------------------------------------------------------------------
      
      if (abs(acc_evid(0) - acc_evid(1)) >= stop_rule) {
        stop("The relative difference in the prior on responses should be smaller than the relative stopping rule before the sampling process. Please adjust \"delta\" or \"prior_on_resp\".");
      }
      
      // for the relative stopping rule, the trial_stim(i) should be either 1 or 2
      if (trial_stim(i) == 1) {
        distr_params = List::create(-1 * discrim/2, distr_add_params(i));
      } else if (trial_stim(i) == 2) {
        distr_params = List::create(discrim/2, distr_add_params(i));
      }
      
      // set the start point
      if (all(is_na(provided_start_point))){ // if users did not provide any start points
        if (i == 0){
          start_point = rDistr(distr_name, distr_params);
          start_point_m = double_to_matrix(start_point, n_chains); // convert start_point to a matrix
          first_sample_idx = 0;
        } else {
          first_sample_idx = 1;
        }
      } else {
        start_point = provided_start_point(i);
        start_point_m = double_to_matrix(start_point, n_chains);
        first_sample_idx = 0;
      }
      
      sampler_results = sampler_mc3_rltv_stop_cpp(
        start_point_m, // NumericMatrix start
        n_chains, // int nChains,
        sigma_prop, // NumericMatrix sigma_prop,
        4, // double delta_T,
        true, // bool swap_all,
        mc3_iterations, // double iterations,
        distr_name, // StringVector distr_name,
        distr_params, // List distr_params,
        false, // bool discreteValues,
        false, // bool isMix,
        1, // NumericVector weights,
        f, // Function custom_func,
        false, // bool useCustom,
        stop_rule, // int stop_rule, // these three parameters are for the ABS
        first_sample_idx, // int save_first_sample
        acc_evid, // NumericVector acc_evid,
        dec_bdry // double dec_bdry
      );
      
      trialSamples = sampler_results["samples"];
      start_point_m = mc3_last_sample(sampler_results["chain"], trialSamples.size()+first_sample_idx-1, n_chains, mc3_iterations);
      trialResponse = sampler_results["response"];
      acc_evid = sampler_results["acc_evid"];
      
      
      trialNDTime = R::runif(nd_time, nd_time+s_nd_time);
      trialDecisionTime = R::rgamma(trialSamples.size(), 1/lambda);
      
      if (trialResponse == 1) {
        conf = acc_evid(0) / sum(acc_evid);
      } else if (trialResponse == 2) {
        conf = acc_evid(1) / sum(acc_evid);
      } else {
        conf = R_NaN;
      }
      
      // update the prior on responses
      acc_evid = clone(prior_on_resp); // reset
      if (prior_depend){
        acc_evid(trial_stim(i)-1) = acc_evid(trial_stim(i)-1) + 1; // update the prior on responses based on the stimulus  
      }
      
      sim_all_trials(i) = List::create(
        _["trial"] = i + 1,
        _["samples"] = trialSamples,
        _["response"] = trialResponse,
        _["stimulus"] = trial_stim(i),
        _["accuracy"] = bool_to_int(trialResponse == trial_stim(i)),
        _["rt"] = trialDecisionTime + trialNDTime,
        _["confidence"] = conf
      );
      break;
    }
  }
  
  return(sim_all_trials);
}
