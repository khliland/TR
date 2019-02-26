#' @title TRYhatCV
#'
#' @description 
#' Function for calculating training set estimates/residuals and Loo estimate/residuals associated
#' with a regularized least squares problem.
#' 
#' @details 
#' Given a TRSVD object, response(s) and regularization parameter value(s), the function finds
#' the fitted Y-values, the LooCV Y-values, the fitted residuals, and the LooCV cross-validated residuals
#' solving the least squares problem implied by the input.
#' 
#' @importFrom pracma repmat
#' 
#' @param decomp TRSVD object created from TRSVDDecomp
#' @param lambdas Vector of candidate regularization parameter values
#' @param Y response vector
#' @param H n*nlambdas matrix of leverage values. Calculated on the fly if not given as input.
#' 
#' @return Returns a list consisting of the fitted responses, the cross-validated responses,
#' the fitted residuals, and the cross-validated residuals.
#' 
#' @seealso 
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}}, 
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#' 
#' @examples 
#' data("gasoline")
#' X        <- as.matrix(gasoline[, 2])
#' Y        <- as.matrix(gasoline[, 1])
#' decomp   <- TRSVDDecomp(X,dtype=1)
#' lambdas  <- 10^seq(-2,5,len=1000)
#' CVOutput <- TRYhatCV(decomp, lambdas, Y)
#' 
#' @export
TRYhatCV <- function(decomp, lambdas, Y, H=TRLeverage(decomp,lambdas)) {
 
  Y        <- as.matrix(Y)
  g        <- dim(Y)[2]
  mY       <- apply(Y,2,mean)
  nlambdas <- length(lambdas)
  n        <- decomp$n
  res      <- array(0,dim=c(n,nlambdas,g))
  rescv    <- array(0,dim=c(n,nlambdas,g))
  Yhat     <- array(0,dim=c(n,nlambdas,g))
  Yhatcv   <- array(0,dim=c(n,nlambdas,g))
  
  Denom <- repmat(lambdas, length(decomp$s), 1)
  Denom <- Denom / t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '/')
  Denom <- Denom + t(repmat(decomp$s, nlambdas, 1)) # sweep(Denom, 1, decomp$s, '+')

  # Append rows to Y matrix if additional criteria have been appended to the data matrix
  # Here mean is appended because it is subtracted again below (yes, it is a 'little' hacky)
  if(decomp$nAddCrit > 0) { Y <- as.matrix(rbind(Y,repmat(mY,decomp$nAddCrit,1))) }  
    
  for(i in 1:g) {
    Yhat[, , i]   <- (mY[i] + decomp$U %*%  ((decomp$s / Denom) * repmat(crossprod(decomp$U, (Y[,i] - mY[i])), 1, nlambdas)))[1:decomp$n]
    #Yhat[, , i]   <- mY[i] + decomp$U %*%  sweep(decomp$s / Denom, 1, t(decomp$U) %*% (Y[,i] - mY[i]), '*')
    res[, , i]    <- Y[1:decomp$n, i] - Yhat[, , i]
    rescv[, , i]  <- res[, , i] / (1-H)
    Yhatcv[, , i] <- Y[1:decomp$n, i] - rescv[, , i]
  }
  
  ret <- list(Yhat=Yhat,Yhatcv=Yhatcv,res=res,rescv=rescv)
}