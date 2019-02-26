#' @title TR SVD Decomposition
#'
#' @description
#' Find SVD of a data matrix to allow for efficient computation of PRESS, GCV, etc for a regularized least squares problem.
#'
#' @details
#' The function finds the singular value decomposition (SVD) of the given data matrix. This decomposition
#' can be used with the other TR functions to efficiently calculate PRESS, GCV, and so on for regularized least squares problems.
#' Additional criteria on the regression coefficients are appended to the data matrix prior to finding the SVD.
#'
#' @importFrom pracma repmat
#' @import SparseM
#' @importFrom Matrix bandSparse
#'
#' @param X Data matrix.
#' @param dtype -1, 0, or positive integer. Controls regularization matrix. -1=Standardization, 0=L_2 regularization,
#' and positive integer gives a discrete derivative matrix of order dtype.
#' @param addCrit Any additional criteria that should be appended to the X matrix. Defaults to FALSE.
#' @param center Boolean. Center data matrix prior to finding the SVD? Defaults to TRUE
#'
#' @return decomp TRSVD object. Contains the SVD, the regularization matrix, and more.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data("gasoline")
#' X      <- as.matrix(gasoline[, 2])
#' decomp <- TRSVDDecomp(X,dtype=1)
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib TR, .registration = TRUE
#' @export
TRSVDDecomp <- function(X, dtype=0, addCrit=FALSE, center=TRUE) {

  X    <- as.matrix(X)
  dims <- dim(X)
  n    <- dims[1]
  p    <- dims[2]

  if(center) {
    mX <- colMeans(X)
    X  <- X - rep(mX, each=n)
  } else {
    mX <- 0
  }

  if(typeof(addCrit) == "list") {
    nAddCrit <- dim(addCrit)[1]
    X <- rbind(as.matrix(X),as.matrix(addCrit))
  } else{
    nAddCrit <- 0
  }

  X <- as.matrix(X)

  if(dtype!=0) {
    L <- TRLmatrix(X,dtype)
    X <- sparseXinvL(X, L)
  } else {
    L <- as.data.frame(list())
  }

  decomp <- dcSVD(X)

  ret <- list(U = decomp$u, s = c(decomp$d), V = decomp$v, L = L, mX = mX,n = n, p = p, centered = center, nAddCrit = nAddCrit)
  class(ret) <- c('TRSVD')
  ret
}
