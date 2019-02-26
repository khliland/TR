#' Specialized X-Y Plotting
#'
#' @aliases summary.TRSVDModel
#' @param x TRSVDModel object (extracted from x if y is missing).
#' @param y not used.
#' @param xvals the coordinates of lines in the plot (\code{xlabs} from y if missing).
#' @param log logarithmic plot, defaults to "x".
#' @param what plot type, defaults to 'prediction', alternatively 'coefficients'.
#' @param which choice of response coefficients if more than one.
#' @param ... arguments to be passed to \code{plot} or \code{lines}.
#' @param object TRSVDModel object.
#'
#' @return
#'
#' @examples
#' data("gasoline")
#' X       <- as.matrix(gasoline[, 2])
#' Y       <- as.matrix(gasoline[, 1])
#' lambdas <- 10^seq(-2,5,len=1000)
#' TRModel <- TRReg(X, Y, dtype=1, lambdas)
#' summary(TRModel)
#' old.par <- par(mfrow = c(2,2), mar = c(3.5,3.5,2.25,3.25), mgp = c(2.25, 1, 0))
#' plot(TRModel)
#' par(old.par)
#'
#' @export
plot.TRSVDModel <- function(x, which=c(1,2,3,4), caption=c("Raw data","Accuracy","Regression coef.","Reference vs predicted"),
                            xvals=NULL, log = "x", main = "", response = 1,
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
  # plot.TRSVDModel <- function(x,y, xvals, log = "x", what = 'prediction', which = 1, ...){
  # Raw spec. | PRESS + GCV
  # -----------------------
  # Reg.coef. | y vs y-hat
  no.xvals <- ifelse(missing(xvals), TRUE, FALSE)
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  if(show[1L]){ # Raw data
    with(x,{
      if(!no.xvals){
        matplot(xvals, t(X), xlab = "abscissa", ylab = "intensity",
                type = 'l', main = caption[1], ...)
      } else {
        if(length(grep("^[-0-9.]+[^0-9]*$", xlabs)) ==
           length(xlabs)) {
          x <- sub("[^0-9]*$", "", xlabs)
          matplot(x, t(X), xlab = "abscissa", ylab = "intensity",
                  type = 'l', main = caption[1], ...)
        } else {
          matplot(t(X), xlab = "abscissa", ylab = "intensity",
                  type = 'l', main = caption[1], ...)
        }
      }
    })
  }

  if(show[2L]){ # PRESS + GCV curves
    with(x,{
      plot(gcv[,response] ~ lambdas, type = 'l', col=2, lty=1, log=log, ylab = "", axes = FALSE, main = caption[2], ...)
      axis(1); axis(2,col.axis = 2); box()
      mtext(side = 2, line = 2, 'GCV', cex = 0.9, col = 2)
      par(new = TRUE)
      plot(rmsecv[,response] ~ lambdas, col = 1, lty = 2, type = "l", log=log, axes = FALSE, bty = "n", xlab = "", ylab = "", ...)
      axis(side=4, at = pretty(range(rmsecv[,response])))
      mtext(side = 4, line = 2, 'RMSECV', cex=0.9)
      # lines(rmsecv[,response] ~ lambdas, col = 1, lty = 2, ...)
      # legend('topleft', lty = 1:2, col = 2:1, legend=c('GCV','RMSECV'))
    })
  }

  if(show[3L]){ # Coefficients (PRESS + GCV optimal)
    with(x,{
      beta  <- c(b[-1,response,])
      betag <- c(bgcv[-1,response,])
      if(!no.xvals){
        x <- xvals
      } else {
        if(length(grep("^[-0-9.]+[^0-9]*$", xlabs)) ==
           length(xlabs)) {
          x <- sub("[^0-9]*$", "", xlabs)
        } else {
          x <- 1:dim(beta)[1]
        }
      }
      plot(beta ~ x, type = 'l', col=2, lty=1, xlab = "abscissa", ylab = "",
           main = caption[3], axes=FALSE, ...)
      lines(x, betag, col=1, lty=2, ...)
      axis(1); axis(2,col.axis = 2); box()
      mtext(side = 2, line = 2, 'GCV', cex = 0.9, col = 2)
      mtext(side = 4, line = 0.5, 'RMSECV', cex = 0.9, col = 1)
      # legend('topleft', lty = 1:2, col = 2:1, legend=c('GCV','RMSECV'))
    })
  }

  if(show[4L]){ # Y vs Y-hat(CV)
    with(x,{
      y <- Y[,response]
      yPRESS <- Y[,response]+residCV$PRESSres[,response]
      yGCV   <- Y[,response]+residCV$GCVres[,response]
      yRange <- range(c(yPRESS, yGCV, y))
      plot(y, yGCV, ylab='Ycv', xlab='Y', col=2, pch='+', main = caption[4],
           xlim = yRange, ylim = yRange, panel.first = lines(yRange,yRange, lty=2, col='#888888'), ...)
      points(y, yPRESS, col=1, pch='x', ...)
      legend('topleft', pch = c('+','x'), col = 2:1, legend=c('GCV','RMSECV'))
    })
  }
}


#' @export
summary.TRSVDModel <- function(object, ...){
  cat("Tikhonov Regularization\n")
  if(object$dtype == 0){
    cat("Ridge")
  } else {
    if(object$dtype > 0)
      cat(paste(object$dtype, c("st","nd","rd",rep("th", max(1, object$dtype)))[object$dtype], " derivative", sep=""))
    else
      cat("standardized Ridge")
  }
  cat("\n\nOptimal fits:\n")
  with(object,
       if(dim(rmsecvOpt)[2] == 1){
         cat(paste("RMSECV = ", round(rmsecvOpt[3],3), " (lambda = ", rmsecvOpt[1], ")\n", sep=""))
         cat(paste("GCV = ", round(gcvOpt[3],3), " (lambda = ", gcvOpt[1], ")\n", sep=""))
       } else {
         for(i in dim(rmsecvOpt)[2]){
           cat(paste("RMSECV = ", round(rmsecvOpt[3,i],3), " (lambda = ", rmsecvOpt[1,i], ")\n", sep=""))
           cat(paste("GCV = ", round(gcvOpt[3,i],3), " (lambda = ", gcvOpt[1,i], ")\n", sep=""))
         }
       })
}

#' @export
print.TRSVDModel <- function(object, ...){
  cat("Tikhonov Regularization\n")
  if(object$dtype == 0){
    cat("Ridge")
  } else {
    if(object$dtype > 0)
      cat(paste(object$dtype, c("st","nd","rd",rep("th", max(1, object$dtype)))[object$dtype], " derivative", sep=""))
    else
      cat("standardized Ridge")
  }
  cat("\n\nOptimal fits:\n")
  with(object,
       if(dim(rmsecvOpt)[2] == 1){
         cat(paste("RMSECV = ", round(rmsecvOpt[3],3), " (lambda = ", rmsecvOpt[1], ")\n", sep=""))
         cat(paste("GCV = ", round(gcvOpt[3],3), " (lambda = ", gcvOpt[1], ")\n", sep=""))
       } else {
         for(i in dim(rmsecvOpt)[2]){
           cat(paste("RMSECV = ", round(rmsecvOpt[3,i],3), " (lambda = ", rmsecvOpt[1,i], ")\n", sep=""))
           cat(paste("GCV = ", round(gcvOpt[3,i],3), " (lambda = ", gcvOpt[1,i], ")\n", sep=""))
         }
       })
}
