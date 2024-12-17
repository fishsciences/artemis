##' Summarize random effects
##'
##' This function returns a summary of the random effects for an
##' eDNA_model produced by \code{eDNA_lmer}
##' 
##' @title Summarize random effect
##' @param object an object of class eDNA_model
##' @param FUN a posterior summary function, default \code{quantile}
##' @param probs probabilities for the posterior summary function
##' @param ... extra args passed to \code{FUN}
##' @return matrix, with one row per random effect, and one column per
##'     probability
##' @author Matt Espe
##' @export
setMethod("ranef", "eDNA_model_lmer",
          function(object, FUN = quantile, probs = c(0.025, 0.5, 0.975), ...)
          {
              rands = extract(object@stanfit, pars = "rand_betas")$rand_betas
              ans = apply(rands, 2, quantile, probs, ...)
              colnames(ans) = colnames(object@random_x)
              t(ans)
          })
