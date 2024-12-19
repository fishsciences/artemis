
##' A 'loo' method for eDNA_model objects, which is simply a wrapper
##' around \code{loo(extract_log_lik(x))}
##'
##' See ?loo::loo for more information about this calculation
##' @title Approximate leave-one-out cross-validation
##' @param x a eDNA_model object
##' @param ... additional args passed to \code{loo}
##' @return a list with the results of the calculation. 
##' @author Matt Espe
##' @export
loo.eDNA_model = function(x,
                          pars = "log_lik",
                          ...,
                          save_psis = FALSE,
                          cores = getOption("mc.cores", 1))
## Based on LOO documentation https://mc-stan.org/loo/reference/loo.html
{
  stopifnot(length(pars) == 1L)
  LLarray <- x@fit$draws(variables = pars)
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(LLarray,
                 r_eff = r_eff,
                 cores = cores,
                 save_psis = save_psis)
}


##' A 'waic' method for eDNA_model objects, which is simply a wrapper
##' around \code{waic(extract_log_lik(x))}
##'
##' See ?loo::waic for more information about this calculation
##' @title Approximate leave-one-out cross-validation
##' @param x a eDNA_model object
##' @param ... additional args passed to \code{waic}
##' @return a list with the results of the calculation. 
##' @author Matt Espe
##' @export
waic.eDNA_model = function(x, ...) {
    waic(extract_log_lik(x@fit), ...)
}
