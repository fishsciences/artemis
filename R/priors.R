##' Normal prior distribution
##'
##' Parameters for the normal distribution, to be used for setting
##' priors on model estimates.
##' @title Normal prior distribution
##' @param location numeric, the mean of the normal distribution
##' @param scale numeric, the sd of the normal distribution
##' @param autoscale logical, whether the priors should be scaled. See
##'     \code{?rstanarm::priors} for details on the scaling.
##' @return named list
##' @author Matt Espe
##' @export
normal = function(location = 0, scale = 1, autoscale = TRUE)
{
    list(dist = "normal", df = NA, location = location,
         scale = scale, autoscale = autoscale)
}
