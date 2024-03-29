% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm.perf.R
\name{gbm.perf}
\alias{gbm.perf}
\title{GBM performance}
\usage{
gbm.perf(object, plot.it = TRUE, oobag.curve = FALSE, overlay = TRUE, method)
}
\arguments{
\item{object}{A \code{\link{gbm.object}} created from an initial call to
\code{\link{gbm}}.}

\item{plot.it}{An indicator of whether or not to plot the performance
measures. Setting \code{plot.it = TRUE} creates two plots. The first plot
plots \code{object$train.error} (in black) and \code{object$valid.error} 
(in red) versus the iteration number. The scale of the error measurement, 
shown on the left vertical axis, depends on the \code{distribution} 
argument used in the initial call to \code{\link{gbm}}.}

\item{oobag.curve}{Indicates whether to plot the out-of-bag performance
measures in a second plot.}

\item{overlay}{If TRUE and oobag.curve=TRUE then a right y-axis is added to
the training and test error plot and the estimated cumulative improvement 
in the loss function is plotted versus the iteration number.}

\item{method}{Indicate the method used to estimate the optimal number of
boosting iterations. \code{method = "OOB"} computes the out-of-bag estimate
and \code{method = "test"} uses the test (or validation) dataset to compute 
an out-of-sample estimate. \code{method = "cv"} extracts the optimal number 
of iterations using cross-validation if \code{gbm} was called with
\code{cv.folds} > 1.}
}
\value{
\code{gbm.perf} Returns the estimated optimal number of iterations.
  The method of computation depends on the \code{method} argument.
}
\description{
Estimates the optimal number of boosting iterations for a \code{gbm} object
and optionally plots various performance measures
}
\seealso{
\code{\link{gbm}}, \code{\link{gbm.object}}
}
\author{
Greg Ridgeway \email{gregridgeway@gmail.com}
}
\keyword{nonlinear}
\keyword{nonparametric}
\keyword{survival}
\keyword{tree}
