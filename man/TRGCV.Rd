% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TRGCV.R
\name{TRGCV}
\alias{TRGCV}
\title{TRGCV}
\usage{
TRGCV(decomp, lambdas, Y, H = TRLeverage(decomp, lambdas),
  centerResponse = TRUE, TRResOutput = TRRes(decomp, lambdas, Y,
  centerResponse = centerResponse))
}
\arguments{
\item{decomp}{TRSVD object created from TRSVDDecomp}

\item{lambdas}{Vector of candidate regularization parameter values}

\item{Y}{response vector}

\item{H}{n*nlambdas matrix of leverage values. Calculated on the fly if not given as input.}

\item{centerResponse}{Boolean. Center response(s)? Defaults to TRUE.}
}
\value{
The output is a list containing (i) the GCV values, (ii) a 3*g matrix of optimal models (for each response the
GCV-minimal lambda-value is given, the index of the GCV-minimal lambda value, as well as the associated GCV value),
(iii) the leverage values. The leverage values are also included in the output list.
}
\description{
Finds the GCV-statistic for the given regularized least squares problem and the minimum GCV-value.
}
\details{
The function returns the GCV-statistic associated with all the regularization parameter values in the lambda vector.
The regularization parameter values minimising the GCV-statistic
is stored in the gcvOpt matrix in the output list, together with the associated index and RMSECV-value.
Leverage values can be given as input or calculated on the fly.
}
\examples{
data(gasoline)
X         <- as.matrix(gasoline[, 2])
Y         <- as.matrix(gasoline[, 1])
decomp    <- TRSVDDecomp(X,dtype=1)
lambdas   <- 10^seq(-2,5,len=1000)
GCVOutput <- TRGCV(decomp, lambdas, Y)

}
\seealso{
\code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
\code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
\code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
\code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
}
