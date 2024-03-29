#' Reconstruct a GBM's Source Data
#' 
#' Helper function to reconstitute the data for plots and summaries. This
#' function is not intended for the user to call directly.
#' 
#' 
#' @param x a \code{\link{gbm.object}} initially fit using \code{\link{gbm}}
#' @return Returns a data used to fit the gbm in a format that can subsequently
#' be used for plots and summaries
#' @author Harry Southworth
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}}
#' @keywords manip
#' @export
reconstructGBMdata <- function(x)
{
   if(!inherits(x, "gbm"))
   {
      stop( "This function is for use only with objects having class 'gbm'" )
   } else
   if (is.null(x$data))
   {
      stop("Cannot reconstruct data from gbm object. gbm() was called with keep.data=FALSE")
   } else
   if (x$distribution$name=="multinomial")
   {
      y <- matrix(x$data$y, ncol=x$num.classes, byrow=FALSE)
      yn <- apply(y, 1, function(z,nc) {(1:nc)[z == 1]},
                  nc = x$num.classes)
      y <- factor(yn, labels=x$classes)
      xdat <- matrix(x$data$x, ncol=ncol(x$data$x.order), byrow=FALSE)
      d <- data.frame(y, xdat)
      names(d) <- c(x$response.name, x$var.names)
   } else
   if (x$distribution$name == "coxph")
   {
      xdat <- matrix(x$data$x, ncol=ncol(x$data$x.order), byrow=FALSE)
      status <- x$data$Misc
      y <- x$data$y[order(x$data$i.timeorder)]
      d <- data.frame(y, status, xdat)
      names(d) <- c(x$response.name[-1], colnames(x$data$x.order))
   }
   else
   {
      y <- x$data$y
      xdat <- matrix(x$data$x, ncol=ncol(x$data$x.order), byrow=FALSE)
      d <- data.frame(y, xdat)
      rn <- ifelse(length(x$response.name) > 1, x$response.name[2], x$response.name)
      names(d) <- c(rn, colnames(x$data$x.order))
   }
   invisible(d)
}
