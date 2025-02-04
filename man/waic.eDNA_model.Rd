% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loo.R
\name{waic.eDNA_model}
\alias{waic.eDNA_model}
\title{Approximate leave-one-out cross-validation}
\usage{
\method{waic}{eDNA_model}(x, pars = "log_lik", ...)
}
\arguments{
\item{x}{a eDNA_model object}

\item{...}{additional args passed to \code{waic}}
}
\value{
a list with the results of the calculation.
}
\description{
A 'waic' method for eDNA_model objects.
}
\details{
See ?loo::waic for more information about this calculation
}
\author{
Matt Espe
}
