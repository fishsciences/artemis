#' The 'artemis' package.
#' 
#' @description The artemis package implements a framework for
#'     simulating data and fitting models for eDNA studies. In
#'     particular, it implements a Bayesian latent-variable model in
#'     which the predictors affect a latent variable, eDNA
#'     concentration. This latent variable is related to an observed
#'     variable, CQ cycles via a standard curve
#'     calibration. CQ cycles, are often truncated at a certain limit,
#'     where more than X cycles is considered a
#'     non-detection.
#'
#' Additional details on the model implemented by the artemis package
#' can be found in the "Getting Started" vignette.
#' 
#' @docType package
#' @name artemis-package
#' @aliases artemis
#' @useDynLib artemis, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import ggplot2
#' @importFrom stats model.matrix model.response model.frame pnorm quantile rnorm terms
#' @importFrom graphics lines segments
#' @importFrom rstan sampling extract
#' @importFrom lme4 lFormula
#' @importFrom rstan sampling
#' @importFrom gridExtra marrangeGrob
#' 
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#' 
NULL
