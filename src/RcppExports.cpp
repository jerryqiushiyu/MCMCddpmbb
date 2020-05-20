// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MCMCddpmbb
Rcpp::List MCMCddpmbb(arma::cube X, arma::Mat<int> dataStructure, double lambda1_value, double lambda0_value, arma::Col<double> alpha_start, arma::Col<double> tau, double W, double m0, double H0, double r0, double s0, double V_start, int sams_iter, int mcmc, int burnin, int thin, int verbose, int seed);
RcppExport SEXP _MCMCddpmbb_MCMCddpmbb(SEXP XSEXP, SEXP dataStructureSEXP, SEXP lambda1_valueSEXP, SEXP lambda0_valueSEXP, SEXP alpha_startSEXP, SEXP tauSEXP, SEXP WSEXP, SEXP m0SEXP, SEXP H0SEXP, SEXP r0SEXP, SEXP s0SEXP, SEXP V_startSEXP, SEXP sams_iterSEXP, SEXP mcmcSEXP, SEXP burninSEXP, SEXP thinSEXP, SEXP verboseSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type dataStructure(dataStructureSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1_value(lambda1_valueSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0_value(lambda0_valueSEXP);
    Rcpp::traits::input_parameter< arma::Col<double> >::type alpha_start(alpha_startSEXP);
    Rcpp::traits::input_parameter< arma::Col<double> >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< double >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< double >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< double >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type V_start(V_startSEXP);
    Rcpp::traits::input_parameter< int >::type sams_iter(sams_iterSEXP);
    Rcpp::traits::input_parameter< int >::type mcmc(mcmcSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMCddpmbb(X, dataStructure, lambda1_value, lambda0_value, alpha_start, tau, W, m0, H0, r0, s0, V_start, sams_iter, mcmc, burnin, thin, verbose, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MCMCddpmbb_MCMCddpmbb", (DL_FUNC) &_MCMCddpmbb_MCMCddpmbb, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_MCMCddpmbb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
