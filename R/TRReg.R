#' @title Tikhonov Regularization Regression
#' @aliases TR
#'
#' @description
#' Function for solving solving regularized least squares problem and finding PRESS and GCV-minimal models.
#'
#' @details
#' The function solves the Tikhonov regularization problem and does model selection by leave-one-out cross-validation
#' (looCV) and generalized cross-validation.
#'
#' @importFrom pracma strcmp
#'
#' @param X Data frame or matrix
#' @param Y Data frame or matrix
#' @param dtype Regularization matrix. -1=Standardization. 0=L_2 regularization. Positive integers give discrete derivative regularization.
#' @param lambdas Vector of candidate regularization parameter values.
#' @param decomp Decomposition used. Only svd is currently supported.
#' @param centerX Center the data matrix? Defaults to TRUE.
#' @param centerY Center the response(s)? Defaults to TRUE.
#' @param addCrit Additional criteria to be appended to the data matrix given as a data frame or matrix. Defaults to FALSE.
#'     ret <- list(b=b, rmsecvOpt=PRESS$rmsecvOpt, rmsecv=PRESS$rmsecv, gcv=GCV$gcv, H=H, decomp=decomp, lambdas=lambdas, bgcv=bgcv, gcvOpt=GCV$gcvOpt)

#' @return The function returns a list containing the following (the name of the object is given in parenthesis after the description):
#' PRESS-minimal regression coefficients for each response (b),
#' a 3*g matrix containing the PRESS-minimal regularization parameter value as well as the associated index in the vector of regularization
#' parameter values and the associated RMSECV-value (rmsecvOpt),
#' a nlambdas*g matrix of RMSECV values for each lambda value and each response (rmsecv),
#' a nlambdas*g matrix of GCV values for each lambda value and each response (gcv),
#' a n*nlambdas matrix of leverage values (H),
#' a TRSVD object consisting of the SVD of the data matrix (right multiplied by the inverse of the regularization matrix) as well as the
#' regularization matrix (decomp),
#' the vector of regularization parameter values (lambdas), GCV-minimal regression coefficients or each response (bgcv),
#' a 3*g matrix containing the GCV-minimal regularization parameter value as well as the associated index in the vector of regularization
#' parameter values and the associated GCV-value (gcvOpt).
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
#'
#' @export
TRReg <- function(X, Y, dtype=0, lambdas=10^seq(-5,20,len=1000), decomp='svd', centerX=TRUE, centerY=TRUE, addCrit=FALSE) {

  if(strcmp(decomp,'svd')){

    decomp   <- TRSVDDecomp(X, dtype, addCrit, centerX)
    H        <- TRLeverage(decomp, lambdas)
    resdenom <- TRRes(decomp, lambdas, Y, centerY)
    PRESS    <- TRPRESS(decomp, lambdas, Y, H, centerY, resdenom)
    GCV      <- TRGCV(decomp, lambdas, Y, H, centerY, resdenom)

    Y    <- as.matrix(Y)
    g    <- dim(Y)[2]

    if(g == 1) {
      b    <- TRbcoefs(decomp, PRESS$rmsecvOpt[1], Y, centerResponse=centerY)
      bgcv <- TRbcoefs(decomp, GCV$gcvOpt[1],      Y, centerResponse=centerY)
    } else {

      b    <- matrix(0, nrow = decomp$p + (centerX | centerY), ncol = g)
      bgcv <- matrix(0, nrow = decomp$p + (centerX | centerY), ncol = g)

      for(i in 1:g) {
        b[, i]    <- TRbcoefs(decomp, PRESS$rmsecvOpt[1,i ], Y[, i], centerResponse=centerY)
        bgcv[, i] <- TRbcoefs(decomp, GCV$gcvOpt[1,i] ,      Y[, i], centerResponse=centerY)
      }
    }

    ret <- list(b=b, rmsecvOpt=PRESS$rmsecvOpt, rmsecv=PRESS$rmsecv, gcv=GCV$gcv,
                H=H, decomp=decomp, lambdas=lambdas, residCV=list(PRESSres=PRESS$residCV, GCVres=GCV$residCV),
                bgcv=bgcv, gcvOpt=GCV$gcvOpt, xlabs=colnames(X), dtype=dtype, X=X,Y=Y)
    class(ret) <- c('TRSVDModel')
    ret
  }
}
