#' @title TRLeverage
#'
#' @description 
#' Calculates leverage values given a SVDDecomp object and a vector of regularization parameter values.
#' 
#' @details 
#' The function uses calculates the leverage values (diagonal elements of the projection matrix) for the
#' regularised least squares problem associated with the given SVDDecomp object and the given vector of
#' regularization parameter values (lambdas). These leverage values can be used to efficiently find
#' the PRESS-statistic and GCV using rank 1 updates.
#' 
#' @importFrom pracma repmat
#' 
#' @param decomp TRSVD object created from TRSVDDEcomp.
#' @param lambdas Vector of regularization parameter values
#' 
#' @return H n*nlambdas matrix of leverage values.
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
#' decomp  <- TRSVDDecomp(X,dtype=1)
#' lambdas <- 10^seq(-2,5,len=1000)
#' H       <- TRLeverage(decomp, lambdas)
#' 
#' @export
TRLeverage <- function(decomp, lambdas) {

  nlambdas <- length(lambdas)
  D        <- replicate(nlambdas, decomp$s^2, 2)
  D        <- D + repmat(lambdas,length(decomp$s),1) # sweep(D,2,lambdas, '+')

  H        <- decomp$U^2 %*% (replicate(nlambdas,decomp$s^2,2)/D)
  
  if(decomp$centered){
    H <- H + 1/decomp$n
  }
  
  H <- H[1:decomp$n,]
  H <- as.matrix(H)
}