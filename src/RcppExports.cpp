// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Majorizor_EM
float Majorizor_EM(arma::mat sig_til, arma::mat sig, arma::mat S, float tau);
RcppExport SEXP _CompStratDesign_Majorizor_EM(SEXP sig_tilSEXP, SEXP sigSEXP, SEXP SSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< float >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(Majorizor_EM(sig_til, sig, S, tau));
    return rcpp_result_gen;
END_RCPP
}
// Majorizor_EM_grad
arma::mat Majorizor_EM_grad(arma::mat sig_til, arma::mat sig, arma::mat S, float tau);
RcppExport SEXP _CompStratDesign_Majorizor_EM_grad(SEXP sig_tilSEXP, SEXP sigSEXP, SEXP SSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< float >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(Majorizor_EM_grad(sig_til, sig, S, tau));
    return rcpp_result_gen;
END_RCPP
}
// Majorizor_linear
float Majorizor_linear(arma::mat sig_til, arma::mat sig, arma::mat S);
RcppExport SEXP _CompStratDesign_Majorizor_linear(SEXP sig_tilSEXP, SEXP sigSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(Majorizor_linear(sig_til, sig, S));
    return rcpp_result_gen;
END_RCPP
}
// Majorizor_linear_grad
arma::mat Majorizor_linear_grad(arma::mat sig_til, arma::mat sig, arma::mat S);
RcppExport SEXP _CompStratDesign_Majorizor_linear_grad(SEXP sig_tilSEXP, SEXP sigSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(Majorizor_linear_grad(sig_til, sig, S));
    return rcpp_result_gen;
END_RCPP
}
// c_prox
arma::vec c_prox(arma::mat S, arma::imat G, double lambda, int maxiter, double tol, int verbose);
RcppExport SEXP _CompStratDesign_c_prox(SEXP SSEXP, SEXP GSEXP, SEXP lambdaSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(c_prox(S, G, lambda, maxiter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// c_ggb_psd
arma::vec c_ggb_psd(arma::mat S, arma::imat G, double lambda, float delta, int maxiter, double tol, int verbose);
RcppExport SEXP _CompStratDesign_c_ggb_psd(SEXP SSEXP, SEXP GSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< float >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(c_ggb_psd(S, G, lambda, delta, maxiter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// c_ggb_mm_linear
arma::vec c_ggb_mm_linear(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose);
RcppExport SEXP _CompStratDesign_c_ggb_mm_linear(SEXP SSEXP, SEXP sig_tilSEXP, SEXP GSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP tSEXP, SEXP BSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(c_ggb_mm_linear(S, sig_til, G, lambda, tau, t, B, maxiter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// c_ggb_mm_EM
arma::vec c_ggb_mm_EM(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose);
RcppExport SEXP _CompStratDesign_c_ggb_mm_EM(SEXP SSEXP, SEXP sig_tilSEXP, SEXP GSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP tSEXP, SEXP BSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_til(sig_tilSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type G(GSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(c_ggb_mm_EM(S, sig_til, G, lambda, tau, t, B, maxiter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CompStratDesign_Majorizor_EM", (DL_FUNC) &_CompStratDesign_Majorizor_EM, 4},
    {"_CompStratDesign_Majorizor_EM_grad", (DL_FUNC) &_CompStratDesign_Majorizor_EM_grad, 4},
    {"_CompStratDesign_Majorizor_linear", (DL_FUNC) &_CompStratDesign_Majorizor_linear, 3},
    {"_CompStratDesign_Majorizor_linear_grad", (DL_FUNC) &_CompStratDesign_Majorizor_linear_grad, 3},
    {"_CompStratDesign_c_prox", (DL_FUNC) &_CompStratDesign_c_prox, 6},
    {"_CompStratDesign_c_ggb_psd", (DL_FUNC) &_CompStratDesign_c_ggb_psd, 7},
    {"_CompStratDesign_c_ggb_mm_linear", (DL_FUNC) &_CompStratDesign_c_ggb_mm_linear, 10},
    {"_CompStratDesign_c_ggb_mm_EM", (DL_FUNC) &_CompStratDesign_c_ggb_mm_EM, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_CompStratDesign(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
