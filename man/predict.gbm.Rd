\name{predict.gbm}
\alias{predict.gbm}
\title{ Predict method for GBM Model Fits }
\description{
  Predicted values based on a generalized boosted model object 
}
\usage{
predict.gbm(object, newdata, n.trees, ...)
}
\arguments{
  \item{object}{ Object of class inheriting from "gbm" (\code{\link{gbm.object}}) }
  \item{newdata}{ Data frame in which to predict }
  \item{n.trees}{ Number of trees used in the prediction }
  \item{\dots}{ further arguments passed to or from other methods }
}
\details{
\code{predict.gbm} runs each observation in \code{newdata} through the trees
described by \code{object$trees[1:n.trees]}, the first n.trees in the boosting
sequence. The function returns sum of the predictions from each of the trees.
}
\value{
Returns a vector of predictions on the scale of f(x). For example, for the
Bernoulli loss the returned value is on the log odds scale, poisson loss on the
log scale, and cox proportional hazard is on the log hazard scale.
}
\author{Greg Ridgeway \email{gregr@rand.org}}
\seealso{
\code{\link{gbm}}, \code{\link{gbm.object}}
}

\keyword{ models }
\keyword{ regression }
