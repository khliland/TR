#' @title TR Plegendre
#'
#' @description
#' Creates Legendre polynomials of given degree and size using a QR factorisation.
#'
#' @details
#' The function creates orthonormal Legendre polynomials of length p and degree d.
#' This is done by sampling the function x^i uniformly at p points in the interval
#' [-1,1]. Orthogonal vectors are then obtained by a QR factorisation.
#'
#' @param d non-negative integer. Degree of polynomials
#' @param p positive integer. Size of vectors
#'
#' @return matrix containing the orthonormal Legendre polynomials.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' Polynomials <- TRPlegendre(2,401)
#'
#' @export
TRPlegendre <- function(d,p) {

  x <- seq(from=-1, to=1, length=p)
  t(qr.Q(qr(sapply(0:d, function(i) x^i))))
  # P <- matrix(1, nrow=p, ncol=d+1)
  # x <- seq(from=-1, to=1, length=p)
  #
  # if(d>0) {
  #   for(i in 1:d) {
  #     P[,i+1] <- x^i
  #   }
  # }
  #
  # qrPol <- t(qr.Q(qr(P)))
  # qrPol <- as.data.frame(t(qrPol))
}
