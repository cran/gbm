## ----shrinkViters, echo=FALSE, out.width="100%", fig.cap="Out-of-sample predictive performance by number of iterations and shrinkage. Smaller values of the shrinkage parameter offer improved predictive performance, but with decreasing marginal improvement.", fig.alt="Out-of-sample predictive performance by number of iterations and shrinkage."----
knitr::include_graphics(if (knitr::is_latex_output()) {
  "shrinkage-v-iterations.pdf"
} else {
  "shrinkage-v-iterations.png"
}, auto_pdf = FALSE)

## ----oobperf, echo=FALSE, out.width="100%", fig.cap="Out-of-sample predictive performance of OOB, test set, and CV methods for selecting the optimal number of iterations. The vertical axis plots performance relative to the best. The boxplots indicate relative performance across thirteen real datasets from the UCI repository. See `demo(OOB-reps)`.", fig.alt="Comparison of methods for estimating the optimal number of boosting iterations."----
knitr::include_graphics(if (knitr::is_latex_output()) {
  "oobperf2.pdf"
} else {
  "oobperf2.png"
}, auto_pdf = FALSE)

