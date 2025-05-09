
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


has_intercept = function(X)
{
    as.integer(any(grepl("\\(Intercept\\)", colnames(X))))
}

conc_thresh = function(a, b, thresh = 40)
{
    (thresh - a) / b
}

##' Convenience function to convert from Cq to [eDNA]
##'
##' Convenience function to convert from Cq to [eDNA]
##' 
##' @title Convert CQ value to [eDNA]
##' @param Cq_values numeric vector, value of CQ
##' @param std_curve_alpha the alpha (intercept) value for the standard curve
##' @param std_curve_beta the beta (slope) value for the standard curve
##' @return numeric vector, [eDNA] values
##' @author Matt Espe
##' @export
cq_to_lnconc = function(Cq_values, std_curve_alpha, std_curve_beta)
{
    (Cq_values - std_curve_alpha) / std_curve_beta
}

##' Convenience function for converting values
##'
##' Convenience function for converting values
##' @title Convert [eDNA] to Cq
##' @param x log[eDNA] values
##' @param std_curve_alpha the alpha (intercept) value for the
##'     standard curve
##' @param std_curve_beta the beta (slope) value for the standard
##'     curve
##' @param censor logical, whether Cq values larger than the upper_Cq
##'     should be censored. If TRUE, these any values above the
##'     threshold will be replaced with the upper_Cq value.
##' @param upper_Cq the max Cq value
##' @return vector of Cq values
##' @author Matt Espe
##' @export
lnconc_to_cq = function(x,
                        std_curve_alpha,
                        std_curve_beta,
                        censor = TRUE, upper_Cq = 40)

{
    ans = std_curve_beta * x +
        std_curve_alpha
    
    if(censor)
        ans[ans > upper_Cq] = upper_Cq
    ans
}


zero_sum_constrain = function(y)
## R version of the sum to zero function inside the Stan models with random effects
{

  N = length(y)
    z = rep(0, N + 1)
    sum_w = 0
    for (ii in 1:N) {
      i = N - ii + 1  
      n = i 
      w = y[i] * (sqrt(n * (n + 1)))^-1 
      sum_w = sum_w + w 
      z[i] = z[i] + sum_w      
      z[i + 1] = z[i + 1] - w * n     
    }
    return(z)
 }
