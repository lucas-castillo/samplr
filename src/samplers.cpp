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

