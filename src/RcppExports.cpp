// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sampler_mcmc_cpp
List sampler_mcmc_cpp(NumericVector start, NumericMatrix sigma_prop, int iterations, String distr_name, List distr_params, bool discreteValues);
RcppExport SEXP _samplr_sampler_mcmc_cpp(SEXP startSEXP, SEXP sigma_propSEXP, SEXP iterationsSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP discreteValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_prop(sigma_propSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< String >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type discreteValues(discreteValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_mcmc_cpp(start, sigma_prop, iterations, distr_name, distr_params, discreteValues));
    return rcpp_result_gen;
END_RCPP
}
// sampler_mc3_cpp
List sampler_mc3_cpp(NumericVector start, int nChains, NumericMatrix sigma_prop, double delta_T, bool swap_all, double iterations, String distr_name, List distr_params, bool discreteValues);
RcppExport SEXP _samplr_sampler_mc3_cpp(SEXP startSEXP, SEXP nChainsSEXP, SEXP sigma_propSEXP, SEXP delta_TSEXP, SEXP swap_allSEXP, SEXP iterationsSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP discreteValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_prop(sigma_propSEXP);
    Rcpp::traits::input_parameter< double >::type delta_T(delta_TSEXP);
    Rcpp::traits::input_parameter< bool >::type swap_all(swap_allSEXP);
    Rcpp::traits::input_parameter< double >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< String >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type discreteValues(discreteValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_mc3_cpp(start, nChains, sigma_prop, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues));
    return rcpp_result_gen;
END_RCPP
}
// sampler_hmc_cpp
List sampler_hmc_cpp(NumericVector start, String distr_name, List distr_params, double epsilon, int L, int iterations);
RcppExport SEXP _samplr_sampler_hmc_cpp(SEXP startSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP epsilonSEXP, SEXP LSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< String >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_hmc_cpp(start, distr_name, distr_params, epsilon, L, iterations));
    return rcpp_result_gen;
END_RCPP
}
// sampler_nuts_cpp
List sampler_nuts_cpp(NumericVector start, String distr_name, List distr_params, double epsilon, int iterations, double delta_max);
RcppExport SEXP _samplr_sampler_nuts_cpp(SEXP startSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP epsilonSEXP, SEXP iterationsSEXP, SEXP delta_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< String >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type delta_max(delta_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_nuts_cpp(start, distr_name, distr_params, epsilon, iterations, delta_max));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_samplr_sampler_mcmc_cpp", (DL_FUNC) &_samplr_sampler_mcmc_cpp, 6},
    {"_samplr_sampler_mc3_cpp", (DL_FUNC) &_samplr_sampler_mc3_cpp, 9},
    {"_samplr_sampler_hmc_cpp", (DL_FUNC) &_samplr_sampler_hmc_cpp, 6},
    {"_samplr_sampler_nuts_cpp", (DL_FUNC) &_samplr_sampler_nuts_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_samplr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}