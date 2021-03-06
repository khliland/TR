% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TRSVCV.R
\name{TRSVCV}
\alias{TRSVCV}
\title{TR Segmented Virtual Cross-Validation}
\usage{
TRSVCV(X, Y, blocks, dtype = 0, lambdas = 10^seq(-5, 20, len = 1000))
}
\arguments{
\item{X}{data frame (predictors).}

\item{Y}{data frame (response(s)).}

\item{blocks}{Either positive integer or a vector. The vector is an indexing vector
for the segments. If all samples are ordered by segments and all segments have the same size
then a positive integer indicating the segment size can be used instead.}

\item{dtype}{Regularization matrix.}

\item{lambdas}{Vector of candidate regularization parameter values.}
}
\value{
Returns a list consisting of (i) SVCV optimal regression coefficients, as well as leverage values,
gcv and rmsecv for all regularization parameter values, and the SVD of the transformed data matrix.
}
\description{
Function does segmented virtual cross-validation for a TR problem.
}
\details{
The function does segmented virtual cross-validation (SVCV) for a TR problem.
The function is useful when one has data that consists of highly similar segments,
such as when obtaining multiple spectra of the same sample.
}
\examples{
data       <- data(fishoil)
X          <- fishoil$Raman
Y          <- fishoil$Iodine
replicates <- fishoil$replicates
lambdas    <- 10^seq(4,15,len=1000)
SVCVOut    <- TRSVCV(X, Y, replicates, 0, lambdas)

}
\seealso{
\code{\link{TRReg}}, \code{\link{TRSVDDecomp}}, \code{\link{TRkFoldCV}}, \code{\link{TRLmatrix}},
\code{\link{TR1SE}}, \code{\link{TRchi2Test}}, \code{\link{TRSVCV}}, \code{\link{TRNumOpt}},
\code{\link{TRPlegendre}}, \code{\link{TRLeverage}}, \code{\link{TRYhatCV}}, \code{\link{TRbcoefs}},
\code{\link{TRPRESS}}, \code{\link{TRGCV}}, \code{\link{TRPred}}, \code{\link{TRDiff}}
}
