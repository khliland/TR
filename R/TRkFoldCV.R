#' @title TR k-fold CV
#'
#' @description 
#' Does k-fold cross-validation for regularized least squares problems.
#' 
#' @details 
#' The function does k-fold cross-validation for regularized least squares problems.
#' The desired regularization matrix is indicated as input together with the candidate regularization parameter values.
#' Multivariate problems are treated as a series of univariate problems.
#' The function finds PRESS-minimal and GCV-minimal regression coefficients.
#' 
#' @param X data frame (predictors)
#' @param Y data frame (response(s))
#' @param dtype Regularization type. -1, 0, or positive integer. -1=Standardization, 0=L_2 regularization,
#' and positive integer gives a discrete derivative matrix of order dtype.
#' @param lambdas Vector of candidate regularization parameter values.
#' @param nFolds Positive integer. Number of folds to be used in cross-validation or
#' and indexing vector for the folds. Defaults to 5.
#' 
#' @return Returns a list consisting of (i) The regression coefficients minimising RMSECV,
#' (ii) RMSECV-values for all the regularization parameter values given as input,
#' (iii) a 3*g matrix of optimal models (for each response the
#' press-minimal lambda-value is given, the index of the press-minimal lambda value, as well as the associated RMSECV value),
#' (iv) TRSVD object associated with the full X-matrix
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
#' TRModel <- TRkFoldCV(X, Y, dtype=1, lambdas, nFolds=5)
#' 
#' @export
TRkFoldCV <- function(X, Y, dtype=0, lambdas=10^seq(-5,20,len=1000), nFolds=5) {
  
  Y        <- as.matrix(Y)
  n        <- nrow(X)
  p        <- ncol(X)
  g        <- ncol(Y)
  nlambdas <- length(lambdas)
  
  # Creating cross-validation folds
  if(length(nFolds) == 1) {
    inds <- sample(n)
    folds <- cut(seq(1,n), breaks=nFolds, labels=FALSE)
  } else {
    folds  <- nFolds
    nFolds <- unique(nFolds)
    inds   <- 1:n
  }
  
  X <- X[inds, ]
  Y <- as.matrix(Y[inds, ])
  
  rmsecv    <- matrix(0, nrow=nlambdas, ncol=g)
  rmsecvOpt <- matrix(0, nrow=3, ncol=g)
  
  for(i in 1:nFolds) {
    testInds <- which(folds == i)
    Xtest    <- X[testInds, ]
    Ytest    <- as.matrix(Y[testInds, ])
    Xtrain   <- X[-testInds, ]
    Ytrain   <- as.matrix(Y[-testInds, ])
    
    decomp <- TRSVDDecomp(Xtrain, dtype)
    bcoefs <- TRbcoefs(decomp, lambdas, Ytrain)
    pred   <- TRPred(Xtest, bcoefs, Ytest)
    rmsecv <- rmsecv + pred[[3]]
  }
  
  decomp <- TRSVDDecomp(X, dtype)
  bcoefs <- matrix(0, nrow=(p+1), ncol=g)
  
  for(i in 1:g) {
    rmsecv[, i]     <- sqrt(1/n * rmsecv[, i])
    rmsecvOpt[2, i] <- which.min(rmsecv[, i])
    rmsecvOpt[1, i] <- lambdas[rmsecvOpt[2, i]]
    rmsecvOpt[3, i] <- rmsecv[rmsecvOpt[2, i], i]
    bcoefs[, i]     <- TRbcoefs(decomp,lambdas=rmsecvOpt[1, i], Y[, i])
  }
  
  ret <- list(bcoefs=bcoefs, rmsecv=rmsecv, rmsecvOpt=rmsecvOpt, decomp=decomp)
}