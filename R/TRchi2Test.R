#' @title TRchi2Test
#'
#' @description
#' Uses a Chi^2-test to select a simpler model than the PRESS-minimal model for regularized least squares problems.
#'
#' @details
#' A Chi^2-test is used to select a simpler model than the PRESS-minimal model for regularized least squares problems
#' by selecting the model with the largest regularization parameter value that is not significantly different
#' (as determined by the significance level) from the PRESS-optimal model.
#'
#' @param decomp Output from TRSVDDecomp
#' @param lambdas Vector of candidate regularization parameter values
#' @param rmsecv RMSECV values corresponding to the lambda vector
#' @param Y Response matrix/vector
#' @param alpha Significance level. Defaults to 0.2.
#'
#' @return rmsecvoptchi2 3*g matrix containing for each response the selected regularization parameter value, the index of this regularization
#' parameter value, as well as the associated RMSECV value.
#'
#' @seealso
#' \code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
#' \code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
#' \code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
#' \code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
#'
#' @examples
#' data("gasoline")
#' X          <- as.matrix(gasoline[, 2])
#' Y          <- as.matrix(gasoline[, 1])
#' lambdas    <- 10^seq(-2,5,len=1000)
#' TRModel    <- TRReg(X, Y, dtype=1, lambdas)
#' chi2Output <- TRchi2Test(TRModel$decomp, lambdas, TRModel$rmsecv, Y, alpha=0.2)
#' bchi2      <- TRbcoefs(TRModel$decomp, chi2Output[1], Y)
#'
#' @export
TRchi2Test <- function(decomp, lambdas, rmsecv, Y, alpha=0.2) {

  Y <-as.matrix(Y)
  g <- dim(Y)[2]

  rmsecvOptchi2 <- matrix(0,3,g)
  nlambdas      <- length(lambdas)
  press         <- rmsecv^2 # * decomp$n

  for(i in 1:g) {
    rmsecvOptchi2[2,i] <- which.min(press[,i])
    rmsecvOptchi2[1,i] <- lambdas[rmsecvOptchi2[2,i]]

    min0 <- min(press[, i]) * decomp$n / qchisq(alpha, decomp$n)
    S    <- (sum(Y[, i]) - Y[, i]) / (decomp$n-1) - Y[, i]
    pressNoModel <- crossprod(S)

    if(pressNoModel>min0) {
      index <- which(press[, i] < min0)
      if(isempty(index)) {

      } else {
        rmsecvOptchi2[2, i] <- index[length(index)]
        rmsecvOptchi2[1, i] <- lambdas[rmsecvOptchi2[2, i]]
      }
    } else {
      rmsecvOptchi2[2, i] <- nlambdas
      rmsecvOptchi2[1, i] <- lambdas[rmsecvOptchi2[2, i]]
    }

    rmsecvOptchi2[3, i] <- rmsecv[rmsecvOptchi2[2, i]]
  }

  rmsecvOptchi2
}
