\name{summary.gbm}
\alias{summary.gbm}
\title{ Summary of a gbm object }
\description{
Computes the relative influence of each variable in the gbm object.
}
\usage{
summary.gbm(object,
            cBars=object$cCols,
            n.trees=object$n.trees,
            plotit=TRUE,
            order=TRUE,
            ...)
}
\arguments{
\item{object}{a \code{gbm} object created from an initial call to 
\code{\link{gbm}}.}  
\item{cBars}{ the number of bars to plot. If \code{order=TRUE} the only the 
variables with the \code{cBars} largest relative influence will appear in the 
barplot. If \code{order=FALSE} then the first \code{cBars} variables will
appear in the plot. In either case, the function will return the relative 
influence of all of the variables.}
\item{n.trees}{ the number of trees used to generate the plot. Only the first
\code{n.trees} trees will be used.}
\item{plotit}{ an indicator as to whether the plot is generated. }
\item{order}{ an indicator as to whether the plotted and/or returned relative 
influences are sorted. }
\item{...}{ other arguments passed to the plot function. }
}
\details{
For \code{distribution="gaussian"} this returns exactly a set of 
Type III sum of squares for each variable normalized to sum to 100. For other
loss functions this returns the reduction attributeable to each varaible in sum 
of squared error in predicting the gradient on each iteration. It describes the 
relative influence of each variable in reducing the loss function. See the 
references below for exact details on the computation.
}
\value{
Returns a data frame where the first component is the variable name and the 
second is the computed relative influence, normalized to sum to 100.
}
\references{ 
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.
}
\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{gbm}} }

\keyword{ hplot }
