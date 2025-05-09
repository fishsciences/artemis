% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_detect.R
\name{est_p_detect_lm}
\alias{est_p_detect_lm}
\alias{est_p_detect_count}
\title{Estimate the probability of detection}
\usage{
est_p_detect_lm(
  variable_levels,
  betas,
  ln_eDNA_sd,
  std_curve_alpha,
  std_curve_beta,
  n_rep = 1:12,
  prob_zero = 0.08,
  model_fit = NULL,
  upper_Cq = 40
)

est_p_detect_count(variable_levels, betas, n_rep = 1:12, model_fit = NULL)
}
\arguments{
\item{variable_levels}{numeric vector, with each element
corresponding to the condition to estimate the probability of
detection.}

\item{betas}{numeric vector, the effect sizes for each of the
variable level}

\item{ln_eDNA_sd}{the measurement error on ln[eDNA]. If a
model_fit is provided and this is missing, the estimated
sd(ln_eDNA) from the model will be used.}

\item{std_curve_alpha}{the alpha for the std. curve formula for
conversion between log(concentration) and CQ}

\item{std_curve_beta}{the alpha for the std. curve formula for
conversion between log(concentration) and CQ}

\item{n_rep}{the number of replicate measurements at the levels
specified}

\item{prob_zero}{the probability of seeing a non-detection,
i.e. zero, from a zero-inflated process. Defaults to 8%, which
is the rate of inflated zeros in a large sampling experiment.}

\item{model_fit}{optional, a model fit from \code{eDNA_lm} or
\code{eDNA_lmer}.  If this is provided, an estimate derived
from the posterior estimates of beta is calculated.}

\item{upper_Cq}{the upper limit on detection. Converted to the
lower_bound of detection internally}
}
\value{
object of class "eDNA_p_detect" with the estimates of the
    probability of detection for the variable levels provided.
}
\description{
Estimate the probability of detection
}
\details{
This function estimates the probability of getting a positive
detection for an eDNA survey given a set of predictors. This can
be useful when trying to take the estimates from a preliminary
study and use those estimates to inform the deployment of future
sampling schemes. The function assumes that you have either an
idea of the effects of the various predictors, for example from a
previous study, or a fit model with estimates of the effect sizes.

This function takes one circumstance at a time, and calculates the
range of outcomes given a number of repeated sampling
attempts. The probability calculated is the probability of getting
at least one positive detection. For details on the underlying
model and assumptions for this calculation, please refer to the
package vignette.
}
\section{Notes on random effects}{


This function deals with random effects in two different
ways. First, when we desire to see the probability of detection
for a specific instance of a random effect, users can specify the
random effect as just another effect by specifying the random
effect = 1 in the variable list, and then the size of the random
effect. However, when users wish to estimate the probability of
detection in cases where random effects are generated from a
distribution of random effects, this can be accomplished by adding
the standard deviation of the random effect to the
\code{Cq_sd}. This takes advantage of the fact that random effects
are just another source of variation, and that sum of random
normal distributions is itself a random normal distribution.
}

\examples{

est_p_detect_lm(variable_levels = c(Intercept = 1, Distance = 100, Volume = 20),
             betas = c(Intercept = -10.5, Distance = -0.05, Volume = 0.001),
             ln_eDNA_sd = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5,
             n_rep = 1:12)

}
\author{
Matt Espe
}
