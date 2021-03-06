% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TR1SE.R
\name{TR1SE}
\alias{TR1SE}
\title{TR 1 standard error}
\usage{
TR1SE(decomp, lambdas, Y)
}
\arguments{
\item{data("gasoline")}{X       <- as.matrix(gasoline[, 2])
Y       <- as.matrix(gasoline[, 1])
decomp  <- TRSVDDecomp(X,dtype=1)}
}
\value{
rmsecvopt1SE 3*g matrix containing for each response the selected regularization parameter value, the index of this regularization
parameter value, as well as the associated RMSECV value.
}
\description{
Uses the "1 standard error rule" to select a simpler model than the PRESS-minimal model for regularized least squares problems.
}
\details{
The function uses the "1 standard error rule" as described  in Elements of Statistical Learning for model selection.
The idea is to treat the squared cross-validation error as a sample and calculate the standard error of this sample.
The "1 standard error rule" then selects the simples model (largest regularization parameter value) that is within one standard
error of the model minimising the PRESS-statistic.
}
\examples{
data("gasoline")
X           <- as.matrix(gasoline[, 2])
Y           <- as.matrix(gasoline[, 1])
decomp      <- TRSVDDecomp(X,dtype=0)
lambdas     <- 10^seq(-4,5,len=1000)
TR1SEOutput <- TR1SE(decomp, lambdas, Y)
b1SE        <- TRbcoefs(decomp, TR1SEOutput[1], Y)

}
\seealso{
\code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
\code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
\code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
\code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
}
