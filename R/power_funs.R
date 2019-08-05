
est_p_detect = function(variable_levels,
                        betas, 
                        Cq_sd,
                        std_curve_alpha, std_curve_beta,
                        n_rep = 1:12,
                        model_fit = NULL,
                        target_detection = 0.8,
                        upper_Cq = 40)

{
    if(!is.null(dim(variablelevels)))
        stop("Sorry, only one set of variable levels at a time currently supported")
    if(is.null(model_fit) && length(variablelevels) != length(betas))
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
