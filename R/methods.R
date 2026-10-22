##' Summarize random effects
##'
##' This function returns a summary of the random effects for a model
##' fit by \code{eDNA_lmer}, \code{eDNA_zinf_lmer}, or
##' \code{eDNA_count_lmer}. Each random effect is summarized by
##' applying \code{FUN} to its posterior draws. Rows are named
##' \code{term:group[level]}, e.g. \code{(Intercept):FilterID[3]}.
##'
##' The estimated standard deviation of each group of random effects
##' is reported by \code{summary}.
##'
##' @title Summarize random effect
##' @param object a model fit with random effects
##' @param FUN a posterior summary function, default \code{quantile}
##' @param probs probabilities for \code{quantile}; ignored for other
##'     functions
##' @param ... extra args passed to \code{FUN}
##' @return matrix, with one row per random effect, and one column per
##'     value returned by \code{FUN}
##' @author Matt Espe
##' @export
setMethod("ranef", "eDNA_model_mixed",
          function(object, FUN = quantile, probs = c(0.025, 0.5, 0.975), ...)
          {
            rands = unclass(object@fit$draws("rand_betas", format = "draws_matrix"))
            ans = if(identical(FUN, quantile))
                      apply(rands, 2, FUN, probs = probs, ...)
                  else apply(rands, 2, FUN, ...)
            ans = if(is.matrix(ans)) t(ans) else as.matrix(ans)
            rownames(ans) = object@random_names
            ans
          })
