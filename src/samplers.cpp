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
