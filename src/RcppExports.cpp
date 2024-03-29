// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sampler_mh_cpp
List sampler_mh_cpp(NumericVector start, NumericMatrix sigma_prop, int iterations, StringVector distr_name, List distr_params, bool discreteValues, bool isMix, NumericVector weights, Function custom_func, bool useCustom);
RcppExport SEXP _samplr_sampler_mh_cpp(SEXP startSEXP, SEXP sigma_propSEXP, SEXP iterationsSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP discreteValuesSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP custom_funcSEXP, SEXP useCustomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_prop(sigma_propSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type discreteValues(discreteValuesSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Function >::type custom_func(custom_funcSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustom(useCustomSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_mh_cpp(start, sigma_prop, iterations, distr_name, distr_params, discreteValues, isMix, weights, custom_func, useCustom));
    return rcpp_result_gen;
END_RCPP
}
// sampler_mc3_cpp
List sampler_mc3_cpp(NumericMatrix start, int nChains, NumericMatrix sigma_prop, double delta_T, bool swap_all, double iterations, StringVector distr_name, List distr_params, bool discreteValues, bool isMix, NumericVector weights, Function custom_func, bool useCustom);
RcppExport SEXP _samplr_sampler_mc3_cpp(SEXP startSEXP, SEXP nChainsSEXP, SEXP sigma_propSEXP, SEXP delta_TSEXP, SEXP swap_allSEXP, SEXP iterationsSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP discreteValuesSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP custom_funcSEXP, SEXP useCustomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_prop(sigma_propSEXP);
    Rcpp::traits::input_parameter< double >::type delta_T(delta_TSEXP);
    Rcpp::traits::input_parameter< bool >::type swap_all(swap_allSEXP);
    Rcpp::traits::input_parameter< double >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type discreteValues(discreteValuesSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Function >::type custom_func(custom_funcSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustom(useCustomSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_mc3_cpp(start, nChains, sigma_prop, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues, isMix, weights, custom_func, useCustom));
    return rcpp_result_gen;
END_RCPP
}
// sampler_hmc_cpp
List sampler_hmc_cpp(NumericVector start, StringVector distr_name, List distr_params, double epsilon, int L, int iterations, bool isMix, NumericVector weights, Function custom_func, bool useCustom);
RcppExport SEXP _samplr_sampler_hmc_cpp(SEXP startSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP epsilonSEXP, SEXP LSEXP, SEXP iterationsSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP custom_funcSEXP, SEXP useCustomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Function >::type custom_func(custom_funcSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustom(useCustomSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_hmc_cpp(start, distr_name, distr_params, epsilon, L, iterations, isMix, weights, custom_func, useCustom));
    return rcpp_result_gen;
END_RCPP
}
// sampler_mc_rec_cpp
List sampler_mc_rec_cpp(NumericVector start, int nChains, double delta_T, bool swap_all, double iterations, StringVector distr_name, List distr_params, bool discreteValues, bool isMix, NumericVector weights, Function custom_func, bool useCustom, double epsilon, int L, double alpha);
RcppExport SEXP _samplr_sampler_mc_rec_cpp(SEXP startSEXP, SEXP nChainsSEXP, SEXP delta_TSEXP, SEXP swap_allSEXP, SEXP iterationsSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP discreteValuesSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP custom_funcSEXP, SEXP useCustomSEXP, SEXP epsilonSEXP, SEXP LSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< double >::type delta_T(delta_TSEXP);
    Rcpp::traits::input_parameter< bool >::type swap_all(swap_allSEXP);
    Rcpp::traits::input_parameter< double >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type discreteValues(discreteValuesSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Function >::type custom_func(custom_funcSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustom(useCustomSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_mc_rec_cpp(start, nChains, delta_T, swap_all, iterations, distr_name, distr_params, discreteValues, isMix, weights, custom_func, useCustom, epsilon, L, alpha));
    return rcpp_result_gen;
END_RCPP
}
// sampler_nuts_cpp
List sampler_nuts_cpp(NumericVector start, StringVector distr_name, List distr_params, double epsilon, int iterations, double delta_max, bool isMix, NumericVector weights, Function custom_func, bool useCustom);
RcppExport SEXP _samplr_sampler_nuts_cpp(SEXP startSEXP, SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP epsilonSEXP, SEXP iterationsSEXP, SEXP delta_maxSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP custom_funcSEXP, SEXP useCustomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type delta_max(delta_maxSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Function >::type custom_func(custom_funcSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustom(useCustomSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler_nuts_cpp(start, distr_name, distr_params, epsilon, iterations, delta_max, isMix, weights, custom_func, useCustom));
    return rcpp_result_gen;
END_RCPP
}
// gridDensity_cpp
NumericVector gridDensity_cpp(StringVector distr_name, List distr_params, bool isMix, NumericVector weights, NumericVector xxRange, NumericVector yyRange, int cellsPerRow, Function densityFunc, bool useCustomDensity);
RcppExport SEXP _samplr_gridDensity_cpp(SEXP distr_nameSEXP, SEXP distr_paramsSEXP, SEXP isMixSEXP, SEXP weightsSEXP, SEXP xxRangeSEXP, SEXP yyRangeSEXP, SEXP cellsPerRowSEXP, SEXP densityFuncSEXP, SEXP useCustomDensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type distr_name(distr_nameSEXP);
    Rcpp::traits::input_parameter< List >::type distr_params(distr_paramsSEXP);
    Rcpp::traits::input_parameter< bool >::type isMix(isMixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xxRange(xxRangeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yyRange(yyRangeSEXP);
    Rcpp::traits::input_parameter< int >::type cellsPerRow(cellsPerRowSEXP);
    Rcpp::traits::input_parameter< Function >::type densityFunc(densityFuncSEXP);
    Rcpp::traits::input_parameter< bool >::type useCustomDensity(useCustomDensitySEXP);
    rcpp_result_gen = Rcpp::wrap(gridDensity_cpp(distr_name, distr_params, isMix, weights, xxRange, yyRange, cellsPerRow, densityFunc, useCustomDensity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_samplr_sampler_mh_cpp", (DL_FUNC) &_samplr_sampler_mh_cpp, 10},
    {"_samplr_sampler_mc3_cpp", (DL_FUNC) &_samplr_sampler_mc3_cpp, 13},
    {"_samplr_sampler_hmc_cpp", (DL_FUNC) &_samplr_sampler_hmc_cpp, 10},
    {"_samplr_sampler_mc_rec_cpp", (DL_FUNC) &_samplr_sampler_mc_rec_cpp, 15},
    {"_samplr_sampler_nuts_cpp", (DL_FUNC) &_samplr_sampler_nuts_cpp, 10},
    {"_samplr_gridDensity_cpp", (DL_FUNC) &_samplr_gridDensity_cpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_samplr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
