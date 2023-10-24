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


// UTILS

NumericVector gradient(dfunc &func, const NumericVector &x, double Temp = 1)
{

  int n_dim = x.size();
  // only one density function so the Jacobian has only one row
  NumericVector returnVector(n_dim);

  double h = 1e-8;
  NumericVector above, below;

  for (int i = 0; i < n_dim; i++){
    above = clone(x);
    below = clone(x);

    above(i) = above(i) + h;
    below(i) = below(i) - h;

    returnVector(i) = ((func(above)/Temp) - (func(below)/Temp)) / (2 * h);
  }
  return returnVector;
}

double dotProduct(const NumericVector &x, const NumericVector &y)
{
  if (x.size() != y.size()){
    stop("Cannot calculate the dot product of vectors of different length");
  }
  double dotProduct = 0;
  for (int i = 0; i < x.size(); i++){
    dotProduct += x(i) * y(i);
  }
  return dotProduct;
}


double joint_d(
    const NumericVector &theta, 
    const NumericVector &momentum, 
    dfunc &log_func, 
    double Temp=1
){
  return (log_func(theta) - (.5 * dotProduct(momentum, momentum) / Temp)) / Temp;
}

NumericVector alpha_trick(
    NumericVector random_jump,
    NumericVector last_jump,
    double alpha
){
  return alpha * last_jump + pow((1 - pow(alpha, 2)), .5) * random_jump;
}



NumericVector autocorrelated_metropolis_step_cpp(
    NumericMatrix &chain, 
    NumericMatrix &proposals, 
    NumericMatrix &jumps, 
    NumericMatrix &true_jumps, 
    const int &currentIndex, 
    const double &last_prob, 
    const NumericMatrix &sigma_prop, 
    dfunc &pdf, 
    const bool &discreteValues, 
    const double &beta, 
    const double &alpha
){
  NumericVector current_x = chain.row(currentIndex - 1);
  NumericVector last_jump = jumps.row(currentIndex - 1);
  NumericVector zeros (current_x.size());
  arma::mat random_jump_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
  NumericVector random_jump = NumericVector(random_jump_.begin(), random_jump_.end());
  
  NumericVector next_jump = alpha_trick(random_jump, last_jump, alpha);
  true_jumps.row(currentIndex) = next_jump;
  NumericVector proposal = current_x + next_jump;
  
  // if dist is discrete round the proposal to nearest int
  if (discreteValues){
    for (int i = 0; i < proposal.length(); i++){
      proposal(i) = round(proposal(i));
    }
  }
  
  // update proposals matrix
  proposals.row(currentIndex) = proposal;
  
  // calculate proposal density
  double prob_prop = pdf(proposal);
  
  // proposal is accepted with probability prob_prop / prob_curr
  if (last_prob != 0){
    double ratio = prob_prop / last_prob;
    
    // The beta parameter (temperature), beta <= 1, 
    // increases the value of the ratio making hotter chains 
    // more likely to accept proposals
    if ((ratio >= 1) || (R::runif(0,1) < pow(ratio, beta))){
      chain.row(currentIndex) = proposal;
      jumps.row(currentIndex) = next_jump;
      return NumericVector::create(prob_prop, 1);
      
    }
  } else if (prob_prop > 0) {
    chain.row(currentIndex) = proposal;
    jumps.row(currentIndex) = next_jump;
    return NumericVector::create(prob_prop, 1);
  }
  chain.row(currentIndex) = current_x;
  jumps.row(currentIndex) = next_jump * -1;
  return NumericVector::create(last_prob, 0);
}


NumericVector propose(
    NumericVector current_x,
    NumericMatrix sigma_prop, 
    double alpha
){
  NumericVector zeros (current_x.size());
  arma::mat perturbance_ = rmvnorm(1, as<arma::vec>(zeros), as<arma::mat>(sigma_prop));
  NumericVector perturbance = NumericVector(perturbance_.begin(), perturbance_.end());
  
  NumericVector proposal = alpha * current_x + pow((1 - pow(alpha, 2)), .5) * perturbance;
  return(proposal);
}


void leapfrog_step_cpp(
    NumericVector &theta, 
    NumericVector &momentum, 
    const double &epsilon, 
    dfunc &log_pdf, 
    const int &L, 
    double Temp=1
)
{
  // start with half step for momentum
  momentum = momentum + (epsilon/2) * gradient(log_pdf, theta, Temp);
  
  // alternate full steps for position and momentum
  for (int i = 0; i<L;i++){
    theta = theta + epsilon * momentum;
    if (i != (L-1)){
      // full step for momentum except for the end of trajectory
      momentum = momentum + epsilon * gradient(log_pdf, theta, Temp);
    }
  }

  // make a half step (instead of a full one) for the momentum at the end
  momentum = momentum + (epsilon/2) * gradient(log_pdf, theta, Temp);

  // negate momentum to make the proposal symmetric

  momentum = -1 * momentum;
}

NumericVector RecycledMomentumUpdate(
    NumericVector momentum, 
    double alpha, 
    double Temp=1
){
  return alpha * momentum + pow((1 - pow(alpha, 2)), .5) * Rcpp::rnorm(momentum.size(), 0, Temp);
}


///'@export
//[[Rcpp::export]]
List sampler_mh_cpp(
    NumericVector start,
    NumericMatrix sigma_prop,
    int iterations,
    StringVector distr_name,
    List distr_params,
    bool discreteValues,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom,
    double alpha=0
)
{
  // Initialize variables ---------------------------------
  LogicalVector acceptances(iterations);
  int n_dim = start.size();
  dfunc pdf = managePDF(distr_name, distr_params, isMix, weights, false, custom_func, useCustom);
  
  NumericMatrix chain(iterations, n_dim);
  NumericMatrix proposals(iterations, n_dim);
  NumericMatrix jumps(iterations, n_dim);
  NumericMatrix true_jumps(iterations, n_dim);
  NumericMatrix ps(1, iterations);
  
  // first row is start
  chain.row(0) = start;
  ps(0,0) = pdf(start);
  
  
  // Run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){
    NumericVector accept;
    if (i == 1){
      accept = autocorrelated_metropolis_step_cpp(
        chain,      // NumericMatrix &chain, 
        proposals, // NumericMatrix &proposals, 
        jumps, // NumericMatrix &jumps, 
        true_jumps, // NumericMatrix &true_jumps, 
        i,        // const int &currentIndex, 
        ps(0,i-1), // const double &lastP, 
        sigma_prop, // const NumericMatrix &sigma_prop, 
        pdf, // dfunc &pdf, 
        discreteValues, // const bool &discreteValues, 
        1, // const double &beta, 
        0 // const double &alpha
      );
    } else{
      accept = autocorrelated_metropolis_step_cpp(
        chain,      // NumericMatrix &chain, 
        proposals, // NumericMatrix &proposals, 
        jumps, // NumericMatrix &jumps, 
        true_jumps, // NumericMatrix &true_jumps, 
        i,        // const int &currentIndex, 
        ps(0,i-1), // const double &lastP, 
        sigma_prop, // const NumericMatrix &sigma_prop, 
        pdf, // dfunc &pdf, 
        discreteValues, // const bool &discreteValues, 
        1, // const double &beta, 
        alpha // const double &alpha
      );
    }
    
    ps(0,i) = accept(0);
    acceptances(i) = (bool)(accept(1));
  }
  
  return List::create(chain, proposals, acceptances, ps, jumps, true_jumps);
}

///'@export
// [[Rcpp::export]]
List sampler_mc3_cpp(
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
    double alpha=0
)
{
  NumericVector acceptances(nChains);
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
  
  // run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){
    
    for (int ch = 0; ch < nChains; ch++){
      NumericVector accept;
      if (i ==1){
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
  }
  
  
  return List::create(
    _["chain"] = chain,
    _["proposals"] = proposals,
    _["beta"] = beta,
    _["swaps"] = swaps,
    _["swap_accepts_attempts"] = NumericVector::create(swap_accepts, swap_attempts),
    _["acceptances"] = acceptances,
    _["jumps"] = jumpsUse,
    _["jumpsHistory"] = jumpsHistory
  );
}




NumericVector drawMomentum(int dimensions){
  NumericVector retV(dimensions);
  for (int i = 0; i < dimensions; i++){
    retV(i) = R::rnorm(0,1);
  }
  return retV;
  
}
///'@export
// [[Rcpp::export]]
List sampler_hmc_cpp(
    NumericVector start,
    StringVector distr_name,
    List distr_params,
    double epsilon,
    int L,
    int iterations,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom
)
{
  // init vars
  dfunc log_pdf = managePDF(distr_name, distr_params, isMix, weights, true, custom_func, useCustom);
  int dim = start.size();
  NumericMatrix chain(iterations, dim);
  NumericMatrix momentums(iterations, dim);
  int acceptances = 0;
  chain.row(0) = start;
  NumericVector momentum;
  
  for (int i = 1; i < iterations; i++){
    // draw a sample of momentum
    
    momentum = drawMomentum(dim);
    
    // initialize vars
    NumericVector theta_prime  = chain.row(i-1);
    NumericVector momentum_prime = clone(momentum);
    
    // leapfrog for each L step
    leapfrog_step_cpp(theta_prime, momentum_prime, epsilon, log_pdf, L);
    
    // Metropolis - Hastings Acceptance, using the joint density of position + momentum
    double top =  exp(joint_d(theta_prime, momentum_prime, log_pdf));
    double bottom =  exp(joint_d(chain.row(i-1), momentum, log_pdf));
    
    double alpha = top/bottom;
    
    if (R::runif(0,1) <= alpha){
      chain.row(i) = theta_prime;
      momentums.row(i) = momentum_prime;
      acceptances++;
    } else{
      chain.row(i) = chain.row(i-1);
      momentums.row(i) = momentums.row(i-1);
    }
  }
  
  return List::create(
    chain, 
    momentums, 
    ((double)(acceptances)/(double)(iterations))
  );
  
}


///'@export
// [[Rcpp::export]]
List sampler_mc_rec_cpp(
    NumericVector start,
    int nChains, // for non-tempered versions of HMC and HOR, set this to one
    double delta_T, // for all chains having the same temp, set to 0
    bool swap_all,
    double iterations,
    StringVector distr_name,
    List distr_params,
    bool discreteValues,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom,
    double epsilon,
    int L,
    double alpha // for HMC, set this to zero
)
{
  if (alpha < -1 || alpha > 1){
    stop("Alpha must be between -1 and +1 (inclusive)");
  }

  // init vars
  NumericVector acceptances(nChains); // counts accept per chain
  int swap_attempts = 0; // counter
  int swap_accepts= 0; // counter
  NumericMatrix swaps(iterations*nChains, 3); //swap History. Columns: iteration, chain1, chain2
  int n_dim = start.length();

  NumericMatrix chain(iterations * nChains, n_dim);
  NumericMatrix proposals(iterations*nChains, n_dim);
  NumericVector beta(nChains); // vector of temperatures
  NumericMatrix ps(nChains, iterations); // densities
  NumericMatrix momentumsHistory(iterations*nChains, n_dim);
  NumericMatrix momentumsUse(iterations*nChains, n_dim);
  NumericVector momentum;
  
  dfunc log_pdf = managePDF(distr_name, distr_params, isMix, weights, true, custom_func, useCustom);
  // init the chain information
  for (int i = 0; i < nChains; i++){
    chain.row(0 + iterations * i) = start;
    beta(i) = (1 + delta_T * i);
  }
  
  int nSwaps; // how many swaps per iteration?
  if (swap_all){
    nSwaps = floor(nChains/2);
  } else {
    nSwaps = 1;
  }
  
  // v holds the chain ids so that they are swapped
  NumericVector v(nChains);
  v = seq(0, nChains - 1);
  
  // run the sampler ------------------------------------------------
  for (int i = 1; i < iterations; i++){
    
    if (i % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    for (int ch = 0; ch < nChains; ch++){
      if (i == 1){
        momentum = drawMomentum(n_dim) * beta(ch);
      } else{
        momentum = momentumsUse.row((i + iterations * ch) - 1);
        momentum = RecycledMomentumUpdate(momentum, alpha, beta(ch));
      }
      
      NumericVector theta_prime  = chain.row((i + iterations * ch) - 1);
      NumericVector momentum_prime = clone(momentum);
      
      // leapfrog for each L step
      leapfrog_step_cpp(theta_prime, momentum_prime, epsilon, log_pdf, L, beta(ch));
      
      // Metropolis - Hastings Acceptance, using the joint density of position + momentum
      
      double top =  exp(joint_d(theta_prime, momentum_prime, log_pdf, beta(ch)));
      double bottom =  exp(joint_d(chain.row((i + iterations * ch) - 1), momentum, log_pdf, beta(ch)));
      
      double alpha = top/bottom;
      if (R::runif(0,1) <= alpha){
        chain.row(i + iterations * ch) = theta_prime;
        momentumsUse.row(i + iterations * ch) = -1 * momentum_prime;
        momentumsHistory.row(i + iterations * ch) = -1 * momentum_prime;
        acceptances(ch) = acceptances(ch) + 1;
      } else{
        chain.row(i + iterations * ch) = chain.row((i + iterations * ch) -1);
        momentumsUse.row(i + iterations * ch) = -1 * momentumsUse.row((i + iterations * ch) -1);
        momentumsHistory.row(i + iterations * ch) = -1 * momentumsUse.row((i + iterations * ch) -1);
      }
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

        double m_pdf = exp(joint_d(chain.row(i + iterations * m), momentumsUse.row(i + iterations * m), log_pdf, beta(m)));
        double n_pdf = exp(joint_d(chain.row(i + iterations * n), momentumsUse.row(i + iterations * n), log_pdf, beta(n)));

        double m_pdf_swapped = exp(joint_d(chain.row(i + iterations * m), momentumsUse.row(i + iterations * m), log_pdf, beta(n)));
        double n_pdf_swapped = exp(joint_d(chain.row(i + iterations * n), momentumsUse.row(i + iterations * n), log_pdf, beta(m)));        // double n_pdf_swapped = exp(joint_d(chain.row(i + iterations * n), momentums.row(i + iterations * n), log_pdf, beta(m)));
        
        double top = m_pdf_swapped * n_pdf_swapped; // pow(m_pdf, beta(n)) * pow(n_pdf, beta(m));
        double bottom = m_pdf * n_pdf; //pow(m_pdf,beta(m)) * pow(n_pdf, beta(n));
        
        
        if ((bottom != 0 && R::runif(0,1) <= top/bottom) || (bottom == 0 && top > 0)){
          
          NumericVector temp = chain.row(i + iterations * m);
          chain.row(i + iterations * m) = chain.row(i + iterations * n);
          chain.row(i + iterations * n) = temp;
          
          NumericVector tempMomentum_m = momentumsUse.row(i + iterations * m);
          NumericVector tempMomentum_n = momentumsUse.row(i + iterations * n);
          
          if (i != (iterations-1)){
            momentumsUse.row(i + iterations * m) = drawMomentum(n_dim) * beta(m);
            momentumsUse.row(i + iterations * n) = drawMomentum(n_dim) * beta(n);
          }
          
          swaps.row(swap_accepts) = NumericVector::create(i+1, m, n);
          swap_accepts++;
        }
      }
    }
  }
  
  
  return List::create(
    _["chain"] = chain,
    _["proposals"] = proposals,
    _["beta"] = beta,
    _["swaps"] = swaps,
    _["swap_accepts_attempts"] = NumericVector::create(swap_accepts, swap_attempts),
    _["acceptances"] = acceptances,
    _["momentums"] = momentumsUse,
    _["momentumsHistory"] = momentumsHistory
  );
}

// // // NUTS // // //

double sampleDirection(){
  if (R::runif(0,1) < .5){
    return -1;
  } else{
    return 1;
  }
}


NumericMatrix build_tree(const NumericVector &theta, const NumericVector &momentum, const double &u, const double  &v, const double &j, const double &epsilon, const double &delta_max, dfunc &log_pdf){
  NumericMatrix tree(7, theta.size());
  // NumericVector temp(theta.size());
  // Rcout << "Starting tree, depth = " << j <<"\n";
  if (j == 0)
  {
    // base case -> Step in direction v
    NumericVector theta_prime = clone(theta);
    NumericVector momentum_prime = clone(momentum);
    // is the new point in the slice?
    leapfrog_step_cpp(theta_prime, momentum_prime, v * epsilon, log_pdf, 1);

    double n_prime = (double)(u <= exp(joint_d(theta_prime, momentum_prime, log_pdf)));
    // is the simulation wildly inaccurate?
    double s_prime = (double)((log(u) - delta_max) < joint_d(theta_prime, momentum_prime, log_pdf));


    tree.row(0) = theta_prime;
    tree.row(1) = momentum_prime;
    tree.row(2) = theta_prime;
    tree.row(3) = momentum_prime;
    tree.row(4) = theta_prime;
    tree(5, 0) = n_prime;
    tree(6, 0) = s_prime;

    // NumericVector temp2 = NumericVector::create(
    //   1,
    //   exp(joint_d(theta_prime, momentum_prime, log_pdf) - joint_d(theta0, momentum0, log_pdf))
    // );
    // temp.fill(min(temp2));
    // tree.row(7) = temp;
    // temp.fill(1);
    // tree.row(8) = temp;
    return tree;
  }
  else
  {
    // recursion -> implicitly build the right and left sub-trees
    NumericMatrix bt = build_tree(theta, momentum, u, v, j-1, epsilon, delta_max, log_pdf);
    NumericVector theta_minus = bt.row(0);
    NumericVector momentum_minus = bt.row(1);
    NumericVector theta_plus = bt.row(2);
    NumericVector momentum_plus = bt.row(3);
    NumericVector theta_prime = bt.row(4);
    double n_prime = bt(5,0);
    double s_prime = bt(6,0);
    // temp = bt.row(7);
    // double alpha_prime = temp(0);
    // temp = bt.row(8);
    // double n_alpha_prime = temp(0);

    // was the stopping criteria met in the first subtree?
    if (s_prime == 1)
    {
      NumericVector theta_prime_2;
      double n_prime_2;
      double s_prime_2;
      // double alpha_prime_2;
      // double n_alpha_prime_2;

      if (v == -1)
      {
        bt = build_tree(theta_minus, momentum_minus, u, v, j-1, epsilon, delta_max, log_pdf);
        theta_minus = bt.row(0);
        momentum_minus = bt.row(1);
        theta_prime_2 = bt.row(4);
        n_prime_2 = bt(5,0);
        s_prime_2 = bt(6,0);
        // temp = bt.row(7);
        // alpha_prime_2 = temp(0);
        // temp = bt.row(8);
        // n_alpha_prime_2 = temp(0);
      }
      else
      {
        bt = build_tree(theta_plus, momentum_plus, u, v, j-1, epsilon, delta_max, log_pdf);
        theta_plus = bt.row(2);
        momentum_plus = bt.row(3);
        theta_prime_2 = bt.row(4);
        n_prime_2 = bt(5,0);
        s_prime_2 = bt(6,0);
        // temp = bt.row(7);
        // alpha_prime_2 = temp(0);
        // temp = bt.row(8);
        // n_alpha_prime_2 = temp(0);
      }
      // which subtree to propagate a sample from
      if ((R::runif(0,1)) <= (n_prime_2 / (n_prime + n_prime_2)))
      {
        theta_prime = clone(theta_prime_2);
      }
      // alpha_prime = alpha_prime + alpha_prime_2;
      // n_alpha_prime = n_alpha_prime + n_alpha_prime_2;
      // update stopping criterion
      NumericVector diff = theta_plus - theta_minus;
      s_prime = s_prime_2 * (double)(dotProduct(diff, momentum_minus) >= 0) * (double)(dotProduct(diff, momentum_plus) >= 0);
      // Update number of valid points
      n_prime = n_prime + n_prime_2;
    }


    tree.row(0) = theta_minus;
    tree.row(1) = momentum_minus;
    tree.row(2) = theta_plus;
    tree.row(3) = momentum_plus;
    tree.row(4) = theta_prime;
    tree(5,0) = n_prime;
    tree(6,0) = s_prime;
    // temp.fill(alpha_prime);
    // tree.row(7) = temp;
    // temp.fill(n_alpha_prime);
    // tree.row(8) = temp;
    return tree;
  }
}

///'@export
//[[Rcpp::export]]
List sampler_nuts_cpp(
    NumericVector start,
    StringVector distr_name,
    List distr_params,
    double epsilon,
    int iterations,
    double delta_max,
    bool isMix,
    NumericVector weights,
    Function custom_func,
    bool useCustom
  )
{
  //// init vars
  dfunc log_pdf = managePDF(distr_name, distr_params, isMix, weights, true, custom_func, useCustom);
  int dim = start.size();
  NumericMatrix chain(iterations, dim);
  chain.row(0) = start;
  NumericMatrix momentums(iterations, dim);

  arma::mat identityMatrix(dim, dim, arma::fill::eye);
  arma::vec zeroes(dim, arma::fill::zeros);
  arma::mat momentum_0 = rmvnorm(1,  zeroes, identityMatrix);
  NumericVector momentum0 = NumericVector(momentum_0.begin(), momentum_0.end());

   // run sampler
   for (double i = 1; i < iterations; i++){
     // resample momentum
     momentum_0 = rmvnorm(1,  zeroes, identityMatrix);
     momentum0 = NumericVector(momentum_0.begin(), momentum_0.end());

     // get a slice u
     double u = R::runif(0, exp(joint_d(chain.row(i-1), momentum0, log_pdf)));

     // init variables
     NumericVector theta_minus = chain.row(i-1);
     NumericVector theta_plus = chain.row(i-1);
     NumericVector momentum_minus = clone(momentum0);
     NumericVector momentum_plus = clone(momentum0);
     double j = 0;
     chain.row(i) = chain.row(i-1);
     double n = 1;
     double s = 1;

     NumericVector theta_prime;
     double n_prime;
     double s_prime;
     // double alpha;
     // double n_alpha;
     // NumericVector temp(theta_minus.size());


     while (s==1)
     {
       // random direction

       if (sampleDirection() == -1)
       {
         NumericMatrix bt = build_tree(theta_minus, momentum_minus, u, -1, j, epsilon, delta_max, log_pdf);

         theta_minus = bt.row(0);
         momentum_minus = bt.row(1);
         theta_prime = bt.row(4);
         n_prime = bt(5,0);
         s_prime = bt(6,0);
         // temp = bt.row(7);
         // alpha = temp(0);
         // temp = bt.row(8);
         // n_alpha = temp(0);
       }
       else
       {
         NumericMatrix bt = build_tree(theta_plus, momentum_minus, u, 1, j, epsilon, delta_max, log_pdf);

         theta_plus = bt.row(2);
         momentum_plus = bt.row(3);
         theta_prime = bt.row(4);
         n_prime = bt(5,0);
         s_prime = bt(6,0);
         // temp = bt.row(7);
         // alpha = temp(0);
         // temp = bt.row(8);
         // n_alpha = temp(0);
       }

       if (s_prime == 1){
         if (R::runif(0,1) <= (n_prime / n)){
           chain.row(i) = theta_prime;
         }

       }
       n = n + n_prime;
       NumericVector diff = theta_plus - theta_minus;
       s = s_prime * (double)(dotProduct(diff, momentum_minus)>= 0) * (double)(dotProduct(diff, momentum_plus)>= 0);

       j++;
     }
     // if (i <= iterations_adapt)
     // {
     //   bar_Hs(i) =
     //     (1 - (1 / (i + t_0))) * bar_Hs(i-1) +
     //     (1 / (i + t_0)) * (delta - alpha / n_alpha);
     //   epsilons(i) = exp(mu - (sqrt(i) / gamma) * bar_Hs(i));
     //   bar_epsilons(i) = exp(
     //     pow(i, kappa * -1) * log(epsilons(i)) +
     //     (1 - pow(i, kappa * -1)) * log(bar_epsilons(i-1))
     //   );
     // }
     // else
     // {
     //   epsilons(i) = bar_epsilons(iterations_adapt);
     // }
   }
    return List::create(chain);
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
