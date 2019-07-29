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


}

prep_model = function(mod_list, alpha, beta, Cq_upper = 40)
{
    model_data = list(N = length(mod_list$y),
                      n_vars = ncol(mod_list$x),
                      X = mod_list$x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper)

    if(is.null(mod_list$groups)){
        b = list(model_data, 
                 has_rand = 0L,
                 n_rand_var = 0L,
                 n_rand_total = 0L, 
                 rand_var_shared = double(0),
                 groups = double(0))

    } else { 
        rand_idx = unlist(relevel_rands(mod_list$groups))
        rand_var_shared = get_shared_rand(mod_list$groups)
        
        b = list(has_rand = 1L,
                 n_rand_var = ncol(mod_list$groups),
                 n_rand_total = max(rand_idx),
                 rand_var_shared = rand_var_shared,
                 groups = rand_idx,
                 rand_sigma = as.array(rand_sd))
        
    }
    
    model_data = c(model_data, b)

    return(model_data)
}
