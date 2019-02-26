#' @title TR 1 standard error
#'
#' @description
#' Uses the "1 standard error rule" to select a simpler model than the PRESS-minimal model for regularized least squares problems.
#'
#' @details
#' The function uses the "1 standard error rule" as described  in Elements of Statistical Learning for model selection.
#' The idea is to treat the squared cross-validation error as a sample and calculate the standard error of this sample.
#' The "1 standard error rule" then selects the simples model (largest regularization parameter value) that is within one standard
#' error of the model minimising the PRESS-statistic.
#'
#' @param
#' data("gasoline")
#' X       <- as.matrix(gasoline[, 2])
#' Y       <- as.matrix(gasoline[, 1])
#' decomp  <- TRSVDDecomp(X,dtype=1)
#'
#' @return rmsecvopt1SE 3*g matrix containing for each response the selected regularization parameter value, the index of this regularization
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
#' X           <- as.matrix(gasoline[, 2])
#' Y           <- as.matrix(gasoline[, 1])
#' decomp      <- TRSVDDecomp(X,dtype=0)
#' lambdas     <- 10^seq(-4,5,len=1000)
#' TR1SEOutput <- TR1SE(decomp, lambdas, Y)
#' b1SE        <- TRbcoefs(decomp, TR1SEOutput[1], Y)
#'
#' @export
TR1SE <- function(decomp, lambdas, Y) {

  Y            <- as.matrix(Y)
  g            <- dim(Y)[2]
  rmsecvOpt1SE <- matrix(0,3,g)
  PRESSStuff   <- TRPRESS(decomp,lambdas,Y)
  nlambdas     <- length(lambdas)

  for(i in 1:g) {
    YhatcvCurrent <- TRYhatCV(decomp,PRESSStuff$rmsecvOpt[1,i],Y[, i])$Yhatcv[, 1, 1]
    SE <- sd((Y[, i] - YhatcvCurrent)^2) / sqrt(decomp$n)

    # Start counting at minimum value and add back index offset
    # If all lambdas satisfy criterion then Inf is returned with a warning.
    # This is dealt with in the next line so the warning is suppressed.
    rmsecvOpt1SE[2, i] <- suppressWarnings(PRESSStuff$rmsecvOpt[2, i] - 2 +
       min(which(PRESSStuff$press[PRESSStuff$rmsecvOpt[2, i]:nlambdas, i] > SE + PRESSStuff$press[PRESSStuff$rmsecvOpt[2,i], i])))
    if(rmsecvOpt1SE[2, i] == Inf) {rmsecvOpt1SE[2, i] <- nlambdas}

    rmsecvOpt1SE[1, i] <- lambdas[rmsecvOpt1SE[2, i]]
    rmsecvOpt1SE[3, i] <- PRESSStuff$rmsecv[rmsecvOpt1SE[2, i], i]
  }

  rmsecvOpt1SE
}
