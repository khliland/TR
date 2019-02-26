#' @title TRRes
#'
#' @description
#' Calculates the residuals for regularized least squares problem.
#'
#' @details
#' Given a TRSVD object, a vector of regularization parameter values, and a response matrix/vector the
#' function calculates the residuals of the regularized least squares problem implied by  the TRSVD object
#' for the regularization parameter values contained in the vector lambdas.
#'
#' @importFrom pracma repmat
#'
#' @param decomp TRSVD object created from TRSVDDecomp
#' @param lambdas Vector of candidate regularization parameter values
#' @param Y response vector
#' @param centerResponse Boolean. Center response(s)? Defaults to TRUE.
#'
#' @return The output is a a list with an array of dimensions (n,nlambdas,g) (number of samples, regularization parameter values, and responses),
#' and a n*nlambdas matrix of denominator factors used by many of the functions in the package.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data(gasoline)
#' X           <- as.matrix(gasoline[, 2])
#' Y           <- as.matrix(gasoline[, 1])
#' decomp      <- TRSVDDecomp(X,dtype=1)
#' lambdas     <- 10^seq(-2,5,len=1000)
#' TRResOutput <- TRRes(decomp, lambdas, Y)
#'
#' @export
TRRes <- function(decomp, lambdas, Y, centerResponse=TRUE) {

Y        <- as.matrix(Y)
n        <- nrow(Y)
g        <- ncol(Y)
nlambdas <- length(lambdas)
resid    <- array(0,c(n,nlambdas,g))

if(centerResponse) {
  mY <- as.matrix(apply(Y,2,mean))
  Y  <- Y - tcrossprod(repmat(1,decomp$n,1), mY)
} else {
  mY <- matrix(0,nrow=1,ncol=g)
}

Denom <- repmat(lambdas, length(decomp$s), 1)
Denom <- Denom / t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '/')
Denom <- Denom + t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '+')

for(i in 1:g) {resid[, , i] <- as.matrix((Y[,i] - decomp$U %*% (replicate(length(lambdas),decomp$s * (crossprod(decomp$U, Y[,i])), 2) / Denom))[1:decomp$n,])}

ret <- list(resid=resid,Denom=Denom)
}
