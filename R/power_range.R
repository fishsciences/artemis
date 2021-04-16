##' @rdname est_power_range_lmer
est_power_range_lm = function(formula, variable_list,
                     betas, sigma_Cq,
                     std_curve_alpha, std_curve_beta,
                     type = c("exclude_zero", "accuracy"),
                     rep_range = seq(2, 20, 2),
                     accuracy_level = 0.2, 
                     conf_level = 0.95,
                     n_sim = 200L,
                     probs = conf_to_probs(conf_level),
                     upper_Cq = 40,
                     verbose = FALSE)
{
    ans = lapply(rep_range, function(i){
        variable_list$rep = seq(i)
        
        est_power_lm(formula, 
                     variable_list, 
                     betas, sigma_Cq, 
                     std_curve_alpha, std_curve_beta, type,
                     accuracy_level, conf_level, n_sim,
                     probs, upper_Cq, expand.grid(variable_list), verbose)
    })
    names(ans) = paste("n_rep:", rep_range)
    ans
}

##' Estimate power for range of reps.
##'
##' This function estimates power for an eDNA sampling study for a
##' range of potential reps.
##' 
##' @title Estimate power for a range of replicates
##' @param formula a model formula, e.g. \code{y ~ x1 + x2}. For
##'     \code{sim_eDNA_lmer}, random intercepts can also be provided,
##'     e.g. \code{ ( 1 | rep ) }.
##' @param variable_list a named list, with the levels that each
##'     variable can take.  Please note that the variables listed in
##'     the formula, including the response variable, must be present
##'     in the variable_list or in the X design matrix.  Extra
##'     variables, i.e. variables which do not occur in the formula,
##'     are ignored.
##' @param betas numeric vector, the beta for each variable in the
##'     design matrix
##' @param sigma_Cq numeric, the measurement error on CQ.
##' @param sigma_rand numeric vector, the stdev for the random
##'     effects. There must be one sigma per random effect specified
##' @param std_curve_alpha the alpha value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param type either "exclude_zero" or "accuracy". Exclude_zero give
##'     the classic power estimate, i.e. whether 0 is in the
##'     confidence interval for the estimate ("significant"). Accuracy
##'     measures whether the estimated betas are within some
##'     percentage of the "true" betas used to simulate the data.
##' @param rep_range vector, a set of the number of iterations to
##'     calculate the power for
##' @param accuracy_level numeric, between 0 and 1. The percent of the
##'     true betas for the accuracy estimate.
##' @param conf_level numeric, between 0 and 1, representing the
##'     percent of the confidence interval to calculate. If
##'     \code{probs} is not provided, then the interval is assumed to
##'     be symetric.
##' @param n_sim integer, the number of simulations to conduct in
##'     order to estimate the power.
##' @param probs probabilities for the calculation of the confidence
##'     intervals. By default, a symetric set of lower and upper
##'     probabilities is constructed by \code{conf_to_probs})
##' @param upper_Cq numeric, the upper limit on CQ detection. Any
##'     value of log(concentration) which would result in a value
##'     greater than this limit is instead recorded as the limit.
##' @param verbose logical, when TRUE output from
##'     \code{rstan::sampling} is written to the console.
##' @return list, with each element a vector for the estimated power
##'     for each parameter for the corresponding number of replicates.
##' @author Matt Espe
est_power_range_lmer = function(formula, variable_list,
                                betas, sigma_Cq,
                                sigma_rand,
                                std_curve_alpha, std_curve_beta,
                                type = c("exclude_zero", "accuracy"),
                                rep_range = seq(2, 20, 2),
                                accuracy_level = 0.2, 
                                conf_level = 0.95,
                                n_sim = 200L,
                                probs = conf_to_probs(conf_level),
                                upper_Cq = 40,
                                verbose = FALSE)
{
    ans = lapply(rep_range, function(i){
        variable_list$rep = seq(i)
        
        est_power_lm(formula, 
                     variable_list, 
                     betas, sigma_Cq, 
                     std_curve_alpha, std_curve_beta, type,
                     accuracy_level, conf_level, n_sim,
                     probs, upper_Cq, verbose)
    })
    names(ans) = paste("n_rep:", rep_range)
    ans
}
