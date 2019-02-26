#' @title TRLmatrix
#'
#' @description
#' Function to construct common regularization matrices used in regularized least squares problems.
#'
#' @details
#' The function can construct common regularization matrices used in regularized least squares problems.
#' Standardization as well as discrete derivative operators are supported.
#' Discrete derivatives are augmented with Legendre polynomials to obtain a full rank matrix. The Legendre
#' polynomials are scaled according to the square root of the epsilon argument.
#' For L_2 regularization an empty data frame is returned.
#'
#' @param X Data matrix or number of columns in data matrix (the matrix itself is necessary when dtype=-1)
#' @param dtype Positive integer. -1 gives standardization matrix, 0 gives empty matrix (for L_2 regularization), and a
#' positive integer gives a discrete derivative matrix of order dtype.
#' @param epsilon If dtype>=1 the extra rows appended to obtain full rank will be scaled by sqrt(epsilon). Defaults to epsilon=1e-10.
#' @param appendLegendre Boolean. Append Legendre polynomials to L matrix to obtain a full rank regularization matrix? Defaults to TRUE.
#'
#' @return Regularization matrix for TR problem
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' L <- TRLmatrix(401, 1, 1e-10)
#'
#' @export
TRLmatrix <- function(X, dtype=0, epsilon=1e-10, appendLegendre=TRUE) {

  # Check if matrix or number of features is given as input
  if(length(dim(X)) > 1) {
    p <- ncol(X)
  } else {
    p <- X
  }

  if(dtype==0) {
    L <- list()
  } else if(dtype > 0) {
    L <- TRDiff(p,dtype)
    if(appendLegendre) {
      L[(p-dtype+1):p,] <- as(sqrt(epsilon) * TRPlegendre(dtype-1,p), "dgCMatrix")
    }
  } else if(dtype < 0) {
    L <- as(diag(apply(X,2,sd)), "dgCMatrix")
  }
  L
}
