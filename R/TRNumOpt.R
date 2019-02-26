#' @title TR Numerical Optimisation
#'
#' @description 
#' Function uses fminbnd to find the regularization parameter value minimising PRESS or GCV for regularized least squares problems.
#' The function only supports univariate responses.
#' 
#' @importFrom pracma fminbnd
#' 
#' @param decomp TRSVD object created from TRSVDdecomp
#' @param y Response vector
#' @param criterion PRESS or GCV. Minimise PRESS or GCV? Defaults to PRESS
#' @param centerResponse Center the response variable? Defaults to TRUE.
#' @param minbnd Lower bound in fminbnd. Defaults to zero.
#' @param maxParam Upper bound in fminbnd. Defaults to 1e20.
#' 
#' @return optLambda Output from fminbnd.
#' 
#' @seealso 
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}}, 
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#' 
#' @examples 
#' data("gasoline")
#' X           <- as.matrix(gasoline[, 2])
#' Y           <- as.matrix(gasoline[, 1])
#' decomp      <- TRSVDDecomp(X,dtype=1)
#' optLambda   <- TRNumOpt(decomp, Y, "PRESS")
#' 
#' @export
TRNumOpt <- function(decomp, y, criterion="PRESS", centerResponse=TRUE, minbnd = 0, maxbnd=1e20) {

  y <- as.matrix(y)
  if(centerResponse) {y <- y - mean(y)}

  # Append zero rows to response additional criteria have been appended to the data matrix
  if(decomp$nAddCrit > 0) { y <- as.matrix(rbind(y,matrix(0,nrow=decomp$nAddCrit,ncol=1))) }
    
  pressValue <- function(lambda) {
    H     <- TRLeverage(decomp, lambda)
    Denom <- matrix(lambda, nrow=length(decomp$s), ncol=1)
    Denom <- sweep(Denom, 1, decomp$s, '/')
    Denom <- sweep(Denom, 1, decomp$s, '+')
    resid <- (y - decomp$U %*% (decomp$s * (crossprod(decomp$U,y)) / Denom))[1:decomp$n,]      
    press <- sum((resid / (1-H))^2)
  }
  
  gcvValue <- function(lambda) {
    H     <- TRLeverage(decomp, lambda)
    Denom <- matrix(lambda, nrow=length(decomp$s), ncol=1)
    Denom <- sweep(Denom, 1, decomp$s, '/')
    Denom <- sweep(Denom, 1, decomp$s, '+')
    resid <- (y - decomp$U %*% (decomp$s * (crossprod(decomp$U,y)) / Denom))[1:decomp$n,]      
    gcv   <- sum((resid / (1-mean(H)))^2)
  }
  
  if(strcmp(criterion,"PRESS")) { optLambda <- fminbnd(pressValue, minbnd, maxbnd)
  } else if (strcmp(criterion,"GCV")) { optLambda <- fminbnd(gcvValue, minbnd, maxbnd)
  }
  
  optLambda
}