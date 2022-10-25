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
// sort_dbl_C
NumericVector sort_dbl_C(NumericVector x);
RcppExport SEXP _ibdsim2_sort_dbl_C(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sort_dbl_C(x));
    return rcpp_result_gen;
END_RCPP
}
// sample_12_C
double sample_12_C();
RcppExport SEXP _ibdsim2_sample_12_C() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(sample_12_C());
    return rcpp_result_gen;
END_RCPP
}
// sample_int_C
NumericVector sample_int_C(double n, double size);
RcppExport SEXP _ibdsim2_sample_int_C(SEXP nSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_C(n, size));
    return rcpp_result_gen;
END_RCPP
}
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
// build_allelemat_C
NumericMatrix build_allelemat_C(NumericVector pos, List haplolist);
RcppExport SEXP _ibdsim2_build_allelemat_C(SEXP posSEXP, SEXP haplolistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< List >::type haplolist(haplolistSEXP);
    rcpp_result_gen = Rcpp::wrap(build_allelemat_C(pos, haplolist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ibdsim2_recombine", (DL_FUNC) &_ibdsim2_recombine, 3},
    {"_ibdsim2_sort_dbl_C", (DL_FUNC) &_ibdsim2_sort_dbl_C, 1},
    {"_ibdsim2_sample_12_C", (DL_FUNC) &_ibdsim2_sample_12_C, 0},
    {"_ibdsim2_sample_int_C", (DL_FUNC) &_ibdsim2_sample_int_C, 2},
    {"_ibdsim2_convert_pos_C", (DL_FUNC) &_ibdsim2_convert_pos_C, 4},
    {"_ibdsim2_build_allelemat_C", (DL_FUNC) &_ibdsim2_build_allelemat_C, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ibdsim2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
