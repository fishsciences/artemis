% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\docType{class}
\name{eDNA_simulation-class}
\alias{eDNA_simulation-class}
\alias{eDNA_simulation_lmer-class}
\alias{eDNA_simulation_lm-class}
\title{eDNA simulation results}
\description{
An S4 object holding the results of \code{sim_eDNA_lm*}.
}
\section{Slots}{

\describe{
\item{\code{ln_conc}}{matrix, simulated log(eDNA concentration)}

\item{\code{Cq_star}}{matrx, simulated CQ star}

\item{\code{formula}}{the formula used}

\item{\code{variable_levels}}{the variable levels used}

\item{\code{betas}}{the effect levels used}

\item{\code{x}}{the model matrix used}

\item{\code{std_curve_alpha}}{the alpha for the std. curve conversion formual used}

\item{\code{std_curve_beta}}{the beta for the std. curve conversion formula used}

\item{\code{upper_Cq}}{the upper limit for CQ}

\item{\code{groups}}{the grouping variables used}

\item{\code{random_sd}}{the stdev of the random effects}
}}

