// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
// convert_pos_C
NumericVector convert_pos_C(NumericVector pos, NumericVector mapFrom, NumericVector mapTo, double extValue);
RcppExport SEXP _ibdsim2_convert_pos_C(SEXP posSEXP, SEXP mapFromSEXP, SEXP mapToSEXP, SEXP extValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mapFrom(mapFromSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mapTo(mapToSEXP);
    Rcpp::traits::input_parameter< double >::type extValue(extValueSEXP);
    rcpp_result_gen = Rcpp::wrap(convert_pos_C(pos, mapFrom, mapTo, extValue));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ibdsim2_recombine", (DL_FUNC) &_ibdsim2_recombine, 3},
    {"_ibdsim2_sortC", (DL_FUNC) &_ibdsim2_sortC, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ibdsim2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
