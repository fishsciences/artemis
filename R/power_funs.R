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

est_power = function(formula, vars_list,
                     betas, sigma_Cq,
                     sigma_rep,
                     std_curve_alpha, std_curve_beta,
                     conf_level = 0.95,
                     n_sim = 200L,
                     upper_Cq = 40,
                     X = expand.grid(vars_list),
                     verbose = FALSE)
{
    
    sims = sim_eDNA_lmer(formula, vars_list,
                         betas, sigma_Cq,
                         sigma_rep,
                         std_curve_alpha, std_curve_beta,
                         n_sim, upper_Cq,
                         X, verbose)

    fit = eDNA_lmer(formula,
                    cbind(expand.grid(vars_list), Cq = sims@Cq_star[1,]),
                    std_curve_alpha = std_curve_alpha,
                    std_curve_beta = std_curve_beta,
                    upper_Cq = upper_Cq, n_chain = 4, iters = 500, verbose = TRUE)



}
