\name{summary.gbm}
\alias{summary.gbm}
\title{ Summary of a gbm object }
\description{
Computes the relative influence of each variable in the gbm object.
}
\usage{
summary.gbm(object,
            cBars=min(10,object$cCols),
            n.trees=object$n.trees,
            plotit=TRUE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{object}{a \code{gbm} object created from an initial call to 
\code{\link{gbm}}.}  
\item{cBars}{ the number of bars to plot. Only the variables with 
the \code{cBars} largest relative influence will appear in the barplot.}
\item{n.trees}{ the number of trees used to generate the plot. Only the first
\code{n.trees} trees will be used.}
\item{plotit}{ an indicator as to whether the plot is generated. }
\item{...}{ other arguments passed to the plot function. }
}
\details{
Details to come later
}
\value{
Returns a data frame where the first component is the variable name and the 
second is the computed relative influence, normalized to sum to 100.
}
\references{ 
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(4).
}
\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{gbm}} }

\keyword{ hplot }
