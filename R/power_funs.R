##' Estimate the probability of detection
##'
##' Estimate the probability of detection
##' 
##' @title Estimate the probability of detection
##' @param variable_levels numeric vector, with each element corresponding
##'     to the condition to estimate the probability of detection.
##' @param betas numeric vector, the effect sizes for each of the variable
##'     level
##' @param Cq_sd numeric, the measurement error on CQ
##' @param std_curve_alpha the alpha for the std. curve formula for conversion
##'     between log(concentration) and CQ
##' @param std_curve_beta the alpha for the std. curve formula for conversion
##'     between log(concentration) and CQ
##' @param n_rep the number of replicate measurements at the levels specified
##' @param model_fit optional, a model fit from \code{eDNA_lm} or \code{eDNA_lmer}.
##'     If this is provided, an estimate derived from the posterior estimates of beta
##'     is calculated.
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

    ln_conc_hat = if(is.null(model_fit)) {
                      variable_levels %*% betas
                  } else {
                      apply(model_fit@betas, 1, function(y) variable_levels %*% y)
                  }

    ln_thresh = (upper_Cq - std_curve_alpha) / std_curve_beta

    ans = sapply(n_rep, function(i) 1 - (pnorm(ln_thresh, ln_conc_hat, Cq_sd) ^ i))

    structure(ans,
              variable_levels = variable_levels,
              reps = n_rep,
              class = c("eDNA_p_detect", class(ans)))
}

##' @rdname est_power_lmer
##' @export
est_power_lm = function(formula, variable_list,
                     betas, sigma_Cq,
                     std_curve_alpha, std_curve_beta,
                     type = c("exclude_zero", "accuracy"),
                     accuracy_level = 0.2, 
                     conf_level = 0.95,
                     n_sim = 200L,
                     probs = c(0.025, 0.975),
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
##' Estimate the power for a given eDNA design and effect sizes.
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
##' @param n_sim integer, the number of simulations to conduct in
##'     order to estimate the power.
##' @param probs probabilities for the calculation of the confidence
##'     intervals.
##' @param upper_Cq numeric, the upper limit on CQ detection. Any value of
##'     log(concentration) which would result in a value greater than this limit is
##'     instead recorded as the limit.
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
##' @export
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

