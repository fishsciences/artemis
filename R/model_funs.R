eDNA_lm = function(formula, data, fit_fun = rstan::sampling,
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 1L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lm(formula, data)
       
    # This works because Stan ignores extra input data
    md = prep_model(ml, std_curve_alpha, std_curve_beta,
                    Cq_upper = upper_Cq)

    md$y = ml$y
    fit = run_model(data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    return(fit)
}

eDNA_lmer = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 1L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lmer(formula, data)
    
    md = prep_model(ml, std_curve_alpha, std_curve_beta,
                    Cq_upper = upper_Cq)

    md$y = ml$y
    fit = run_model(data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    return(fit)
}

run_model = function(model = stanmodels$eDNA_lm, data, n_chain,
                  iters, verbose, ...)

{
    if(!verbose) sink(tempfile())
    
    fit = sampling(model, data, chains = n_chain,
                  iter = iters,
                  refresh = ifelse(verbose, 100, -1), show_messages = verbose, ...)

    if(!verbose) sink()
    
    cl = if(is_lme4(formula)) "eDNA_model_lmer" else "eDNA_model_lm"
    
    fit = as(fit, "eDNA_model")
    class(fit) = cl
    return(fit)
}
