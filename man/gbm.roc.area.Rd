% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ir.measures.R
\name{gbm.roc.area}
\alias{gbm.roc.area}
\alias{gbm.conc}
\alias{ir.measure.conc}
\alias{ir.measure.auc}
\alias{ir.measure.mrr}
\alias{ir.measure.map}
\alias{ir.measure.ndcg}
\alias{perf.pairwise}
\title{Compute Information Retrieval measures.}
\usage{
gbm.roc.area(obs, pred)

gbm.conc(x)

ir.measure.conc(y.f, max.rank = 0)

ir.measure.auc(y.f, max.rank = 0)

ir.measure.mrr(y.f, max.rank)

ir.measure.map(y.f, max.rank = 0)

ir.measure.ndcg(y.f, max.rank)

perf.pairwise(y, f, group, metric = "ndcg", w = NULL, max.rank = 0)
}
\arguments{
\item{obs}{Observed value.}

\item{pred}{Predicted value.}

\item{x}{Numeric vector.}

\item{y, y.f, f, w, group, max.rank}{Used internally.}

\item{metric}{What type of performance measure to compute.}
}
\value{
The requested performance measure.
}
\description{
Functions to compute Information Retrieval measures for pairwise loss for a
single group. The function returns the respective metric, or a negative
value if it is undefined for the given group.
}
\details{
For simplicity, we have no special handling for ties; instead, we break ties
randomly. This is slightly inaccurate for individual groups, but should have
only a small effect on the overall measure.

\code{gbm.conc} computes the concordance index: Fraction of all pairs (i,j)
with i<j, x[i] != x[j], such that x[j] < x[i]

If \code{obs} is binary, then \code{gbm.roc.area(obs, pred) =
gbm.conc(obs[order(-pred)])}.

\code{gbm.conc} is more general as it allows non-binary targets, but is
significantly slower.
}
\references{
C. Burges (2010). "From RankNet to LambdaRank to LambdaMART: An
Overview", Microsoft Research Technical Report MSR-TR-2010-82.
}
\seealso{
\code{\link{gbm}}
}
\author{
Stefan Schroedl
}
\keyword{models}
