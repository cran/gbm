\name{calibrate.plot}
\alias{calibrate.plot}
\title{Calibration plot}
\description{
An experimental diagnostic tool that plots the fitted values versus the actual average values.
Currently developed for only \code{distribution="bernoulli"}.
}
\usage{
calibrate.plot(y, p, ...)
}
\arguments{
  \item{y}{ the outcome 0-1 variable }
  \item{p}{ the estimated probabilities }
  \item{...}{ other graphics parameters passed on to the plot function }
}
\details{
Uses \code{\link[modreg]{loess}} to estimate E(y|p). Well-calibrated probabilities
imply that E(y|p) = p. The plot also includes a pointwise 95% confidence
band.
}
\value{
\code{calibrate.plot} returns no values.
}
\references{ 
J.F. Yates (1982). "External correspondence: decomposition of the mean
probability score," Organisational Behaviour and Human Performance 30:132-156.

D.J. Spiegelhalter (1986). "Probabilistic Prediction in Patient Management
and Clinical Trials," Statistics in Medicine 5:421-433.
}
\author{Greg Ridgeway \email{gregr@rand.org}}
\examples{
library(rpart)
data(kyphosis)
y <- as.numeric(kyphosis$Kyphosis)-1
x <- kyphosis$Age
glm1 <- glm(y~poly(x,2),family=binomial)
p <- predict(glm1,type="response")
calibrate.plot(y, p, xlim=c(0,0.6), ylim=c(0,0.6))
}
\keyword{ hplot }
