eDNA_lm = function(formula, data, fit_fun = rstan::sampling,
                   n_chain = 1L, iters = 500L, verbose = FALSE, ...)
{
    
    ml = if(is_lme4(formula)) {
             gen_model_list_lm(formula, X)
         } else {
             gen_model_list_lmer(formula, X)
         }
    
    md = prep_sim(ml, std_curve_alpha, std_curve_beta, sigma_Cq, betas)

    if(!verbose) sink(tempfile())
    
    fit = fit_fun(stanmodels$eDNA_sim, data = md, chains = n_chain,
                  algorithm = "Fixed_param", iter = iters,
                  refresh = ifelse(verbose, 100, -1), show_messages = verbose)

    if(!verbose) sink()
    cl = if(is_lme4) "eDNA_model_lm" else "eDNA_model_lmer"
    
    fit = as(fit, cl)
    return(fit)
}

