// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
// we need R.h to manage RNG when repeated calls to R functions (see customPDF)
#include <R.h>

// 
#include "pdf_manage.h"

using namespace Rcpp;

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

NumericVector drawMomentum(int dimensions){
  NumericVector retV(dimensions);
  for (int i = 0; i < dimensions; i++){
    retV(i) = R::rnorm(0,1);
  }
  return retV;
  
}
