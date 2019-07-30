eDNA_lm = function(formula, data, fit_fun = rstan::sampling,
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 1L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lm(formula, data)
       
    # This works because Stan ignores extra input data
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                    Cq_upper = upper_Cq, type = "model")
    md$y = ml$y
    fit = run_model(data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    fit = load_slots_model(fit)
    class(fit) = "eDNA_model_lm"
    
    return(fit)
}

eDNA_lmer = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 1L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lmer(formula, data)
    
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, type = "model")
    md$y = ml$y

    fit = run_model(data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    fit = load_slots_model(fit)
    class(fit) = "eDNA_model_lmer"
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
    
    fit = as(fit, "eDNA_model")
    return(fit)
}

load_slots_model = function(obj)
    {
    obj@formula = get("formula", parent.frame())
    obj@x = as.data.frame(get("md", parent.frame())$X)
    obj@std_curve_alpha = get("std_curve_alpha", parent.frame())
    obj@std_curve_beta = get("std_curve_beta", parent.frame())
    obj@upper_Cq = get("upper_Cq", parent.frame())
    
    obj
}
