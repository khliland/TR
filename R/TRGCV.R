#' @title TRGCV
#'
#' @description
#' Finds the GCV-statistic for the given regularized least squares problem and the minimum GCV-value.
#'
#' @details
#' The function returns the GCV-statistic associated with all the regularization parameter values in the lambda vector.
#' The regularization parameter values minimising the GCV-statistic
#' is stored in the gcvOpt matrix in the output list, together with the associated index and RMSECV-value.
#' Leverage values can be given as input or calculated on the fly.
#'
#' @importFrom pracma repmat
#'
#' @param decomp TRSVD object created from TRSVDDecomp
#' @param lambdas Vector of candidate regularization parameter values
#' @param Y response vector
#' @param H n*nlambdas matrix of leverage values. Calculated on the fly if not given as input.
#' @param centerResponse Boolean. Center response(s)? Defaults to TRUE.
#'
#' @return The output is a list containing (i) the GCV values, (ii) a 3*g matrix of optimal models (for each response the
#' GCV-minimal lambda-value is given, the index of the GCV-minimal lambda value, as well as the associated GCV value),
#' (iii) the leverage values. The leverage values are also included in the output list.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data(gasoline)
#' X         <- as.matrix(gasoline[, 2])
#' Y         <- as.matrix(gasoline[, 1])
#' decomp    <- TRSVDDecomp(X,dtype=1)
#' lambdas   <- 10^seq(-2,5,len=1000)
#' GCVOutput <- TRGCV(decomp, lambdas, Y)
#'
#' @export
TRGCV <- function(decomp, lambdas, Y, H=TRLeverage(decomp,lambdas), centerResponse=TRUE, TRResOutput=TRRes(decomp, lambdas, Y, centerResponse=centerResponse)) {

  Y <- as.matrix(Y)
  g <- dim(Y)[2]

  if(centerResponse) {
    mY <- as.matrix(apply(Y,2,mean))
    Y  <- Y - tcrossprod(repmat(1,decomp$n,1), mY)
  } else {
    mY <- matrix(0,nrow=1,ncol=g)
  }

  # Append zero rows to Y matrix if additional criteria have been appended to the data matrix
  if(decomp$nAddCrit > 0) { Y <- as.matrix(rbind(Y,matrix(0,nrow=decomp$nAddCrit,ncol=g))) }

  nlambdas <- length(lambdas)
  gcv      <- matrix(nrow=nlambdas, ncol=g)
  gcvOpt   <- matrix(nrow=3, ncol=g)
  residCV  <- matrix(nrow=dim(Y)[1], ncol=g)

  # Denom <- repmat(lambdas, length(decomp$s), 1)
  # Denom <- Denom / t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '/')
  # Denom <- Denom + t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '+')

  Denom <- TRResOutput$Denom
  resid <- TRResOutput$resid

  for(i in 1:g) {
    #resid    <- as.matrix((Y[,i] - decomp$U %*% (replicate(length(lambdas),decomp$s * crossprod(decomp$U, Y[,i]), 2) / Denom))[1:decomp$n,])
    residCVall  <- resid[, , i] / (1 - repmat(apply(H,2,mean),decomp$n, 1))
    gcv[, i] <- apply( (residCVall)^2, 2, sum)
    #gcv[,i] <- sqrt(1/decomp$n * gcv[,i])

    gcvOpt[2,i] <- which.min(gcv[,i])
    gcvOpt[1,i] <- lambdas[gcvOpt[2,i]]
    gcvOpt[3,i] <- gcv[gcvOpt[2,i],i]
    residCV[,i] <- residCVall[,gcvOpt[2, i]]
  }

  ret <- list(gcv=gcv, gcvOpt=gcvOpt, H=H, residCV=residCV)
}
