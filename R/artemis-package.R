#' The 'artemis' package.
#' 
#' @description Simulating data and fitting models for eDNA studies
#'     can be accomplished using the artemis package. In particular,
#'     artemis implements a Bayesian latent-variable model in which the
#'     predictors affect a latent variable, eDNA concentration. This
#'     latent variable is related to an observed variable, CQ cycles
#'     via a standard curve calibration. CQ cycles, are often
#'     truncated at a certain limit, where more than X cycles is
#'     considered a non-detection.
#'
#' Additional details on the model implemented by the artemis package
#' can be found in the "Getting Started" vignette.
#' 
#' @docType package
#' @name artemis-package
#' @aliases artemis
#' #@useDynLib artemis, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import ggplot2
#' @import cmdstanr
#' @importFrom stats model.matrix model.response model.frame pnorm quantile rnorm terms sd
#' @importFrom graphics lines segments plot
#' @importFrom rstan sampling extract plot
#' @importFrom lme4 lFormula ranef
#' @importFrom rstan sampling
#' @importFrom loo loo waic extract_log_lik
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#' 
NULL
