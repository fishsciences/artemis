
##' A 'loo' method for eDNA_model objects.
##'
##' See ?loo::loo for more information about this calculation
##' @title Approximate leave-one-out cross-validation
##' @param x a eDNA_model object
##' @param ... additional args passed to \code{$loo()}
##' @return a list with the results of the calculation. 
##' @author Matt Espe
##' @export
loo.eDNA_model = function(x,
                          ...)
## Based on LOO documentation https://mc-stan.org/loo/reference/loo.html
{
  x@fit$loo(...)
}


##' A 'waic' method for eDNA_model objects.
##'
##' See ?loo::waic for more information about this calculation
##' @title Approximate leave-one-out cross-validation
##' @param x a eDNA_model object
##' @param ... additional args passed to \code{waic}
##' @return a list with the results of the calculation. 
##' @author Matt Espe
##' @export
waic.eDNA_model = function(x, ...) {
  LLarray <- x@fit$draws(variables = pars)
  waic(LLarray, ...)
}
