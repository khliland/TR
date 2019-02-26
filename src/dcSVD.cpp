#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List dcSVD(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, X, "both", "dc");
  return List::create(Rcpp::Named("u") = U,
                      Rcpp::Named("d") = S,
                      Rcpp::Named("v") = V);
}

