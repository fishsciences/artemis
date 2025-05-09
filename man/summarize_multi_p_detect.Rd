% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_detect.R
\name{summarize_multi_p_detect}
\alias{summarize_multi_p_detect}
\title{Summarize multiple probability of detection}
\usage{
summarize_multi_p_detect(x, probs = c(0.025, 0.5, 0.975))
}
\arguments{
\item{x}{list, the results from \code{est_p_detect_lm()} when called
with multiple standard curve parameters}

\item{probs}{vector, probabilities for quantiles passed to
\code{quantile()}}
}
\value{
data.frame
}
\description{
Summarizes estimates where multiple standard curves are used in
\code{est_p_detect_lm()}, which results in a list being returned for
each standard curve parameter.
}
\author{
Matt Espe
}
