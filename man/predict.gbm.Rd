\name{predict.gbm}
\alias{predict.gbm}
\title{ Predict method for GBM Model Fits }
\description{
  Predicted values based on a generalized boosted model object 
}
\usage{
predict.gbm(object,newdata,n.trees,
            type="link",
            single.tree=FALSE,
            ...)
}
\arguments{
  \item{object}{ Object of class inheriting from (\code{\link{gbm.object}}) }
  \item{newdata}{ Data frame of observations for which to make predictions }
  \item{n.trees}{ Number of trees used in the prediction }
  \item{type}{ The scale on which gbm makes the predictions }
  \item{single.tree}{If \code{single.tree=TRUE} then \code{predict.gbm} returns
                     only the predictions from tree \code{n.trees}.}
  \item{\dots}{ further arguments passed to or from other methods }
}
\details{
\code{predict.gbm} runs each observation in \code{newdata} through the trees
described by \code{object$trees[1:n.trees]}, the first n.trees in the boosting
sequence. The function returns sum of the predictions from each of the trees.

If \code{object} was fit using \code{\link{gbm.fit}} there will be no 
\code{Terms} component. Therefore, the user has greater responsibility to make
sure that \code{newdata} is of the same format (order and number of variables)
as the one originally used to fit the model.
}
\value{

Returns a vector of predictions. By default the predictions are on the scale of f(x). For example, for the Bernoulli loss the returned value is on the log odds scale, poisson loss on the log scale, and coxph is on the log hazard scale.

If \code{type="response"} then \code{gbm} converts back to the same scale as the outcome. Currently the only effect this will have is returning probabilities for bernoulli and expected counts for poisson. For the other distributions "response" and "link" return the same.
}
\author{Greg Ridgeway \email{gregr@rand.org}}
\seealso{
\code{\link{gbm}}, \code{\link{gbm.object}}
}

\keyword{ models }
\keyword{ regression }
