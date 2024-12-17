##' Normal prior distribution
##'
##' Parameters for the normal distribution, to be used for setting
##' priors on model estimates.
##'
##' These are styled after the distributions in rstanarm.
##' 
##' @title Prior distributions
##' @param location numeric, the mean of the  distribution
##' @param rate numeric, the rate of the exponential distribution
##' @param scale numeric, the sd/scale of the  distribution
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

##' @rdname normal
##' @export
exponential = function(rate = 1, autoscale = FALSE) 
{
    stopifnot(length(rate) == 1)
    list(dist = "exponential", df = NA, location = NA, scale = 1/rate, 
        autoscale = autoscale)
}
