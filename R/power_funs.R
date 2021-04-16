
##' @rdname est_power_lmer
est_power_lm = function(formula, variable_list,
                     betas, sigma_Cq,
                     std_curve_alpha, std_curve_beta,
                     type = c("exclude_zero", "accuracy"),
                     accuracy_level = 0.2, 
                     conf_level = 0.95,
                     n_sim = 200L,
                     probs = conf_to_probs(conf_level),
                     upper_Cq = 40,
                     X = expand.grid(variable_list),
                     verbose = FALSE)
{
    if(!type %in% c("exclude_zero", "accuracy"))
        stop("Invalid type: ", type)
    
    sims = sim_eDNA_lm(formula, variable_list,
                       betas, sigma_Cq,
                       std_curve_alpha, std_curve_beta,
                       n_sim, upper_Cq,
                       X, verbose)

    ans = lapply(seq(n_sim), function(i) {
        fit = eDNA_lm(formula,
                      cbind(sims@x, Cq = sims@Cq_star[i,]),
                      std_curve_alpha = std_curve_alpha,
                      std_curve_beta = std_curve_beta,
                      upper_Cq = upper_Cq, n_chain = 4, iters = 500)
        as.data.frame(t(apply(fit@betas, 2, quantile, probs)),
                      row.names = colnames(fit@x))
    })

    res = if(type == "exclude_zero"){
              lapply(ans, excludes_zero)
          } else { 
              lapply(ans, function(x) within_accuracy(x, betas, accuracy_level))
          }
    
    res = do.call(cbind, res)
    
    return(rowSums(res)/ncol(res))
}

##' Estimate the power for a given eDNA design
##'
##' These functions estimate the power for a given eDNA design and
##' specific effect sizes. The functions allow the user to specify to
##' either use the classic definition of power (i.e., whether an
##' estimate includes 0 in the confidence/credible interval), or
##' accuracy, defined as the estimates being within a certain
##' percentage of the "true" betas. Both can be useful in different
##' circumstances, but eDNA survey studies produce data where accuracy
##' can be a more useful metric. Since the response variable, Cq
##' values, are truncated at an upper limit, it is often possible to
##' detect a "significant effect" which results in Cq values above
##' this upper limit, but the effects are estimated with a high amount
##' of uncertainty. When are primarily interested in the precision of
##' the estimate, the classic definition of power is not helpful. For
##' these reasons, these functions allow the user to specify
##' "accuracy" as the metric of interest.
##' 
##' @title Estimate power
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
##' @param X optional, a design matrix. By default, this is created
##'     from the variable_list using \code{expand.grid()}, which
##'     creates a balanced design matrix. However, the user can
##'     provide their own \code{X} as well, in which case the
##'     variable_list is ignored. This allows users to provide an
##'     unbalanced design matrix.
##' @param verbose logical, when TRUE output from
##'     \code{rstan::sampling} is written to the console.
##' @return a named vector, with one estimate of the power for each
##'     parameter estimate.
##' @author Matt Espe
est_power_lmer = function(formula, variable_list,
                          betas, sigma_Cq,
                          sigma_rand,
                          std_curve_alpha, std_curve_beta,
                          type = c("exclude_zero", "accuracy"),
                          accuracy_level = 0.2, 
                          conf_level = 0.95,
                          n_sim = 200L,
                          probs = conf_to_probs(conf_level),
                          upper_Cq = 40,
                          X = expand.grid(variable_list),
                          verbose = FALSE)
{
    if(!type %in% c("exclude_zero", "accuracy"))
        stop("Invalid type: ", type)
    
    sims = sim_eDNA_lmer(formula, variable_list,
                         betas, sigma_Cq, sigma_rand,
                         std_curve_alpha, std_curve_beta,
                         n_sim, upper_Cq,
                         X, verbose)

    ans = lapply(seq(n_sim), function(i) {
        fit = eDNA_lmer(formula,
                        cbind(sims@x, Cq = sims@Cq_star[i,]),
                        std_curve_alpha = std_curve_alpha,
                        std_curve_beta = std_curve_beta,
                        upper_Cq = upper_Cq, n_chain = 4, iters = 500)
        as.data.frame(t(apply(fit@betas, 2, quantile, probs)),
                      row.names = colnames(fit@x))
    })
    
    res = if(type == "exclude_zero"){
              lapply(ans, excludes_zero)
          } else { 
              lapply(ans, function(x) within_accuracy(x, betas, accuracy_level))
          }
    
    res = do.call(cbind, res)
    
    return(rowSums(res)/ncol(res))
}

