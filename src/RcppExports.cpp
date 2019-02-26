// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// dcSVD
List dcSVD(const arma::mat& X);
RcppExport SEXP _TR_dcSVD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(dcSVD(X));
    return rcpp_result_gen;
END_RCPP
}
// sparseXinvL
NumericMatrix sparseXinvL(NumericMatrix X, SEXP L);
RcppExport SEXP _TR_sparseXinvL(SEXP XSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseXinvL(X, L));
    return rcpp_result_gen;
END_RCPP
}
// sparseinvLX
NumericMatrix sparseinvLX(SEXP L, NumericMatrix X);
RcppExport SEXP _TR_sparseinvLX(SEXP LSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseinvLX(L, X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TR_dcSVD", (DL_FUNC) &_TR_dcSVD, 1},
    {"_TR_sparseXinvL", (DL_FUNC) &_TR_sparseXinvL, 2},
    {"_TR_sparseinvLX", (DL_FUNC) &_TR_sparseinvLX, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
