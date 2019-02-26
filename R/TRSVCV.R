#' @title TR Segmented Virtual Cross-Validation
#'
#' @description
#' Function does segmented virtual cross-validation for a TR problem.
#'
#' @details
#' The function does segmented virtual cross-validation (SVCV) for a TR problem.
#' The function is useful when one has data that consists of highly similar segments,
#' such as when obtaining multiple spectra of the same sample.
#'
#' @importFrom pracma repmat
#'
#' @param X data frame (predictors).
#' @param Y data frame (response(s)).
#' @param blocks Either positive integer or a vector. The vector is an indexing vector
#' for the segments. If all samples are ordered by segments and all segments have the same size
#' then a positive integer indicating the segment size can be used instead.
#' @param dtype Regularization matrix.
#' @param lambdas Vector of candidate regularization parameter values.
#'
#' @return Returns a list consisting of (i) SVCV optimal regression coefficients, as well as leverage values,
#' gcv and rmsecv for all regularization parameter values, and the SVD of the transformed data matrix.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data       <- data(fishoil)
#' X          <- fishoil$Raman
#' Y          <- fishoil$Iodine
#' replicates <- fishoil$replicates
#' lambdas    <- 10^seq(4,15,len=1000)
#' SVCVOut    <- TRSVCV(X, Y, replicates, 0, lambdas)
#'
#' @export
TRSVCV <- function(X, Y, blocks, dtype=0, lambdas=10^seq(-5,20,len=1000)) {

  blockORTH <- function(X, blocks) {
    # Calculates the orthogonal transformation needed for pre-processing
    n            <- nrow(X)
    uniqueBlocks <- unique(blocks)
    U            <- matrix(0, nrow=n, ncol=n)

    for(k in uniqueBlocks) {
      ind <- (blocks == k)
      svdDecomp <- svd(X[ind, ])
      U[ind, ind] <- svdDecomp$u
    }
    U
  }

  X         <- as.matrix(X)
  Y         <- as.matrix(Y)
  n         <- nrow(X)
  p         <- ncol(X)
  g         <- ncol(Y)
  nlambdas  <- length(lambdas)
  bcoefs    <- matrix(0, nrow=(p+1), ncol=g)
  bcoefsgcv <- matrix(0, nrow=(p+1), ncol=g)

  if(length(blocks) == 1) {
    nblocks <- blocks
    blocks  <- rep(1:(nrow(X)/nblocks), each=nblocks)
  }

  Ublocks <- blockORTH(X, blocks)
  bs      <- as.matrix(apply(Ublocks,2,sum)^2) # Leverage values

  # Centering
  mY <- as.matrix(apply(Y,2,mean))
  Y  <- Y - repmat(1,n,g) %*% t(mY)
  mX <- as.matrix(apply(X,2,mean))
  X  <- X - repmat(1,n,1) %*% t(mX)
  # Note that as the centering is done here no centering is done in the function calls below.

  # Transforming data
  X <- crossprod(Ublocks, X)
  Y <- crossprod(Ublocks, Y)

  # Finding PRESS and GCV
  decomp   <- TRSVDDecomp(X,dtype,center=FALSE)
  H        <- TRLeverage(decomp, lambdas)
  H        <- H + repmat(bs/n, 1, nlambdas) # sweep(H, 1, bs/n)
  PRESSOut <- TRPRESS(decomp, lambdas, Y, H, centerResponse=FALSE)
  GCVOut   <- TRGCV(decomp, lambdas, Y, H, centerResponse=FALSE)

  # Calculating regression coefficients
  for(i in 1:g) {
    bcoefs[2:(p+1), i]    <- as.matrix(TRbcoefs(decomp, PRESSOut$rmsecvOpt[1, i], Y[, i], centerResponse=FALSE))
    bcoefsgcv[2:(p+1), i] <- as.matrix(TRbcoefs(decomp, GCVOut$gcvOpt[1, i], Y[, i], centerResponse=FALSE))

    bcoefs[1, i]    <- mY[, i] - t(mX) %*% as.matrix(bcoefs[2:(p+1), i])
    bcoefsgcv[1, i] <- mY[, i] - t(mX) %*% bcoefsgcv[2:(p+1), i]
  }

  ret <- list(bcoefs=bcoefs, rmsecvOpt=PRESSOut$rmsecvOpt, rmsecv=PRESSOut$rmsecv, gcv=GCVOut$gcv, decomp=decomp, bcoefsgcv=bcoefsgcv, gcvOpt=GCVOut$gcvOpt, H=H)
}
