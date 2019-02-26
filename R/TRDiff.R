#' @title TRDiff
#'
#' @description
#' Function for creating square discrete derivative matrix of given size and degree.
#'
#' @details
#' The function is used to create square discrete derivative matrices to be used as regularization matrices in regularized
#' least squares problems.
#' The function is similar to "diff" in MATLAB.
#'
#' @param p Size of matrix
#' @param dtype Derivative order
#'
#' @return sparse Matrix of discrete derivatives.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' L <- TRDiff(401, 1)
#'
#' @export
TRDiff <- function(p, dtype) {
  dPattern <- rev(diff(c(rep(0,dtype),1,rep(0,dtype)),1,dtype))
  pPattern <- 0:dtype
  bandSparse(p, p, pPattern,
                  diagonals = matrix(dPattern, p,length(dPattern), byrow=TRUE),
                  symmetric = FALSE)
  # L <- diag(p+dtype)
  # L <- diff(L,differences=dtype)[,-c((ncol(L)-dtype+1:(ncol(L))))]
  # L <- as.data.frame(L)
}
