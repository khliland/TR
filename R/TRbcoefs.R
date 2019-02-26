#' @title TRbcoefs
#'
#' @description
#' Calculates regression coefficients solving a regularized least squares problem.
#'
#' @details
#' The function calculates regression coefficients for the regularized least squares problem(s) given by the
#' TRSVD object, the response(s) and the regularization paramter value(s).
#'
#' @importFrom pracma isempty
#' @importFrom Matrix bandSparse
#'
#' @param decomp TRSVD object created from TRSVDDecomp.
#' @param lambdas Vector of candidate regularization parameter values.
#' @param Y response(s).
#' @param centerResponse Boolean. Center response(s)? Defaults to TRUE.
#'
#' @return Returns an array of size (p+1)*nlambdas*g containing the regression coefficients.
#' The first element of the first dimension is the constant term.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data("gasoline")
#' X       <- as.matrix(gasoline[, 2])
#' Y       <- as.matrix(gasoline[, 1])
#' decomp  <- TRSVDDecomp(X,dtype=1)
#' bcoef   <- TRbcoefs(decomp, 10, Y)
#'
#' @export
TRbcoefs <- function(decomp, lambdas, Y, centerResponse=TRUE) {

  Y <- as.matrix(Y)
  g <- dim(Y)[2]

  if(centerResponse) {
    mY <- as.matrix(apply(Y,2,mean))
    Y  <- Y - repmat(1,decomp$n,1) %*% t(mY)
  } else {
    mY <- matrix(0,nrow=1,ncol=dim(Y)[2])
  }

  # Append zero rows to Y matrix if additional criteria have been appended to the data matrix
  if(decomp$nAddCrit > 0) { Y <- as.matrix(rbind(Y,matrix(0,nrow=decomp$nAddCrit,ncol=g))) }

  nlambdas <- length(lambdas)
  g        <- dim(Y)[2]
  p        <- decomp$p
  bcoefs   <- array(0,dim=c(p+1,nlambdas,g))

  Denom <- repmat(lambdas, length(decomp$s), 1)
  Denom <- Denom / t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '/')
  Denom <- Denom + t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '+')

  # Finding regression coefficients (except for the constant term)
  for(i in 1:g) {
    bcoefs[2:(p+1), , i] <- decomp$V %*% (1 / (Denom / repmat(crossprod(decomp$U,Y[,i]), 1, nlambdas)))
    #bcoefs[2:(p+1), , i] <- decomp$V %*% (1/sweep(Denom,1,t(decomp$U) %*% Y[,i], '/'))
  }

  # Transform back if regularization matrix is not the identity
  if(!isempty(decomp$L)) {
    for(i in 1:g) {
      bcoefs[2:(p+1), , i] <- sparseinvLX(decomp$L,as.matrix(bcoefs[-1, , i]))
    }
  }

  # Finding constant term if response or data matrix is centered
  for(i in 1:g) {
    if(decomp$centered) {bcoefs[1, , i] <- crossprod(decomp$mX,bcoefs[2:(p+1), , i])}
    if(centerResponse) {bcoefs[1, , i] <- mY[i] - bcoefs[1, , i]}
  }

  if(!decomp$centered & !centerResponse) {
    bcoefs <- bcoefs[2:(p+1), , i]
  }

  bcoefs
}
