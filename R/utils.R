
excludes_zero = function(ci)
{
    ci[2] < 0 | ci[1] > 0
}

within_accuracy = function(ci, true_betas, acc = 0.20)
{
    ci[1] > true_betas * (1 - acc) &
        ci[2] < true_betas * (1 + acc)
}
##' Confidence level to probability
##'
##' 
##' This is a convience function to convert from a confidence level to
##' the lower and upper probability of the confidence interval.
##' @title Confidence to probability
##' @param conf numeric between 0 and 1
##' @return vector of length 2
##' @author Matt Espe
##' @export
##' @examples
##'
##' conf_to_probs(0.95)
conf_to_probs = function(conf)
{
    x = (1 - conf) / 2
    c(x, 1 - x)
}


in_interval = function(ci, beta)
{
    beta > ci[1] & beta < ci[2]
}

has_intercept = function(x)
{
    as.integer(any(colnames(x) == "(Intercept)"))
}
