#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericMatrix sparseXinvL(NumericMatrix X, SEXP L){
  using Eigen::Map;
  using Eigen::MatrixXd;
  typedef Eigen::MappedSparseMatrix<double> MSpMat;
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::SparseLU<SpMat> SpChol;
  typedef Map<MatrixXd> MapMatd;
  const SpMat Lt(as<MSpMat>(L).adjoint());
  const MatrixXd Xt(as<MapMatd>(X).adjoint());
  const SpChol Ch(Lt);
  if (Ch.info() != Eigen::Success)
    return R_NilValue;
  MatrixXd ret = Ch.solve(Xt).adjoint();
  return wrap(ret);
}

// [[Rcpp::export]]
NumericMatrix sparseinvLX(SEXP L, NumericMatrix X){
  using Eigen::Map;
  using Eigen::MatrixXd;
  typedef Eigen::MappedSparseMatrix<double> MSpMat;
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::SparseLU<SpMat> SpChol;
  typedef Map<MatrixXd> MapMatd;
  const SpMat Lt(as<MSpMat>(L));
  const MatrixXd Xt(as<MapMatd>(X));
  const SpChol Ch(Lt);
  if (Ch.info() != Eigen::Success)
    return R_NilValue;
  MatrixXd ret = Ch.solve(Xt);
  return wrap(ret);
}
