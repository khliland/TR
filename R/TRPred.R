#' @title TR Prediction
#'
#' @description 
#' Function for prediction given predictors and regression coefficients. If the true response
#' is given as input then the prediction error will also be calculated.
#' 
#' @importFrom pracma repmat
#' 
#' @param X data frame (predictors).
#' @param bcoefs Regression coefficients.
#' @param Y (Optional) True response(s). Defaults to FALSE
#' 
#' @return If the true response is not given as input the function returns a list consisting
#' of Yhat calculated from the given data and regression coefficients. If the true response is
#' given as input then the output list also contains the RMSEP and squared prediction error.
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
#' lambdas <- 10^seq(-2,5,len=1000)
#' TRModel <- TRReg(X, Y, dtype=1, lambdas)
#' TRPred  <- TRPred(X, TRModel$b, Y)
#' 
#' @export
TRPred <- function(X, bcoefs, Y=FALSE) {
 
  dimB     <- dim(bcoefs)
  nlambdas <- dimB[2]
  g        <- dimB[3]
  if(is.na(g)) { g <- 1}
  n        <- nrow(X)
  Yhat     <- array(0,dim=c(n,nlambdas,g))
  
  # Check if there is a constant term and separating this from the rest of the regression coefficients
  if(dimB[1]-1 == dim(X)[2]) {
    constTerm <- bcoefs[1, , , drop=FALSE]
    bcoefs    <- bcoefs[2:dimB[1], , , drop=FALSE]
  } else {
    constTerm <- matrix(0,nrow=1,ncol=g)
  }
  
  # Finding fitted values
  for(i in 1:g) {
    Yhat[, , i] <- as.matrix(X) %*% as.matrix(bcoefs[, , i]) + tcrossprod(repmat(1,n,1), as.matrix(constTerm[, , i]))
  }
  
  # Calculating prediction error if true response is given as input
  if(is.list(Y) | is.matrix(Y)) {
    rmsep   <- matrix(0, nrow=nlambdas, ncol=g)
    predErr <- matrix(0, nrow=nlambdas, ncol=g)
    
    for(i in 1:g) {
      squareError  <- (as.matrix(Yhat[, , i]) - repmat(as.matrix(Y[, i]), 1, nlambdas))^2
      predErr[, i] <- apply(squareError, 2, sum)
      rmsep[, i]   <- sqrt(predErr[, i])
    }
    
    ret <- list(Yhat=Yhat, rmsep=rmsep, predErr=predErr)
  } else {
    ret <- list(Yhat=Yhat)
  }
  
}