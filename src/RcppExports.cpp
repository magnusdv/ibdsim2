// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// recombine
NumericMatrix recombine(NumericMatrix strand1, NumericMatrix strand2, NumericVector cross);
RcppExport SEXP _ibdsim2_recombine(SEXP strand1SEXP, SEXP strand2SEXP, SEXP crossSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type strand1(strand1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type strand2(strand2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cross(crossSEXP);
    rcpp_result_gen = Rcpp::wrap(recombine(strand1, strand2, cross));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ibdsim2_recombine", (DL_FUNC) &_ibdsim2_recombine, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ibdsim2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
