# scripts to work with the probability of detection either to work
# backwards to find the variable level which corresponds to a certain
# p(detect) for a single sample or the joint p(detect) for a bunch of
# samples

##' Estimate the probability of detection
##'
##' This function estimates the probability of getting a positive
##' detection for an eDNA survey given a set of predictors. This can
##' be useful when trying to take the estimates from a preliminary
##' study and use those estimates to inform the deployment of future
##' sampling schemes. The function assumes that you have either an
##' idea of the effects of the various predictors, for example from a
##' previous study, or a fit model with estimates of the effect sizes.
##'
##' This function takes one circumstance at a time, and calculates the
##' range of outcomes given a number of repeated sampling
##' attempts. The probability calculated is the probability of getting
##' at least one positive detection. For details on the underlying
##' model and assumptions for this calculation, please refer to the
##' package vignette.
##'
##' @section Notes on random effects:
##'
##' This function deals with random effects in two different
##' ways. First, when we desire to see the probability of detection
##' for a specific instance of a random effect, users can specify the
##' random effect as just another effect by specifying the random
##' effect = 1 in the variable list, and then the size of the random
##' effect. However, when users wish to estimate the probability of
##' detection in cases where random effects are generated from a
##' distribution of random effects, this can be accomplished by adding
##' the standard deviation of the random effect to the
##' \code{Cq_sd}. This takes advantage of the fact that random effects
##' are just another source of variation, and that sum of random
##' normal distributions is itself a random normal distribution.
##'
##' @title Estimate the probability of detection
##' @param variable_levels numeric vector, with each element
##'     corresponding to the condition to estimate the probability of
##'     detection.
##' @param betas numeric vector, the effect sizes for each of the
##'     variable level
##' @param Cq_sd numeric, the measurement error on CQ. If a model_fit
##'     is provided and this is missing, the estimated sd(CQ) from the
##'     model will be used.
##' @param std_curve_alpha the alpha for the std. curve formula for
##'     conversion between log(concentration) and CQ
##' @param std_curve_beta the alpha for the std. curve formula for
##'     conversion between log(concentration) and CQ
##' @param n_rep the number of replicate measurements at the levels
##'     specified
##' @param model_fit optional, a model fit from \code{eDNA_lm} or
##'     \code{eDNA_lmer}.  If this is provided, an estimate derived
##'     from the posterior estimates of beta is calculated.
##' @param upper_Cq the upper limit on detection
##' @return object of class "eDNA_p_detect" with the estimates of the
##'     probability of detection for the variable levels provided.
##' @author Matt Espe
##'
##' @examples
##'
##' est_p_detect(variable_levels = c(Intercept = 1, Distance = 100, Volume = 20),
##'              betas = c(Intercept = -10.5, Distance = -0.05, Volume = 0.001),
##'              Cq_sd = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5,
##'              n_rep = 1:12)
##'
##' @export
est_p_detect = function(variable_levels,
                        betas, 
                        Cq_sd,
                        std_curve_alpha, std_curve_beta,
                        n_rep = 1:12,
                        model_fit = NULL,
                        upper_Cq = 40)

{
    if(!is.null(dim(variable_levels)))
        stop("Sorry, only one set of variable levels at a time currently supported")
    if(is.null(model_fit) && length(variable_levels) != length(betas))
        stop("Variable levels and betas cannot be of different lengths")
    if(missing(betas) && is.null(model_fit))
        stop("Must provide either a set of beta values or a model_fit object")
    
    if(!missing(betas)) {
        if(!is.null(model_fit)) dup_arg_warn("beta")
        ln_conc_hat = variable_levels %*% betas
    } else {
        # model_fit provided
        if(!missing(Cq_sd)){
            dup_arg_warn("Cq_sd")
        } else {
            Cq_sd = model_fit@sigma_Cq
        }
        if(!missing(std_curve_alpha) || !missing(std_curve_beta)) {
            dup_arg_warn("std_curve parameters")
        } else {
            std_curve_alpha = model_fit@std_curve_alpha
            std_curve_beta = model_fit@std_curve_beta            
        }
        ln_conc_hat = apply(model_fit@betas, 1, function(y) variable_levels %*% y)
    }
    
    Cq_hat = ln_conc_hat * std_curve_beta + std_curve_alpha
    ans = sapply(n_rep, function(i) 1 - ((1 - pnorm(upper_Cq, Cq_hat, Cq_sd)) ^ i))

    structure(ans,
              variable_levels = variable_levels,
              reps = n_rep,
              class = c("eDNA_p_detect", class(ans)))
}

dup_arg_warn = function(arg)
{
    warning(sprintf("Both %s and model_fit provided. Using %s provided.", arg, arg))
}

if(FALSE){
    ## not ready yet
var_level_given_detect = function(p, model_fit,
                                  varying = character(), fixed_levels,
                                  upper_Cq = model_fit@upper_Cq,
                                  Cq_sd = model_fit@Cq_sd,
                                  std_curve_alpha = model_fit@std_curve_alpha,
                                  std_curve_beta = model_fit@std_curve_beta)
    # replicated the analytic steps taken in Expt5Design.R
{
    ln_conc_thresh = (upper_Cq - std_curve_alpha) / std_curve_beta

    if(is.null(names(fixed_levels)))
        stop("Please provide a named vector of fixed levels")
    
    fixed_order = match(names(fixed_levels), colnames(model_fit@x))

    if(any(is.na(fixed_order)))
        stop("Unable to match names of 'fixed_levels' to parameters")

    var = match(varying, colnames(model_fit@x))

    if(any(is.na(var)))
        stop("Could not match 'varying' with a parameter. Please check the spelling.")
    
    # Fix this
    var_at_threshold = (ln_conc_thresh -
                        (fixed_levels[fixed_order] %*%
                         model_fit@x[,names(fixed_levels[fixed_order])]) /
                        model_fix@betas)
    
}
}
