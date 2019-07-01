#  Functions to call Stan model and calculate power


sim_experiment = function(vars_list, betas, sigma_rand,
                          sigma_conc, sigma_Cq, ...)
{
}

prep_model_data = function(vars_list, betas, sigma_rand,
                           sigma_Cq,
                           X = expand.grid(vars_list),
                           rand_var_names = c("tech_rep", "rep"))
{
    tmp = sim_data(, betas, sigma_rand, sigma_Cq, X, rand_var_names)
    rand_idx = relevel_rands(tmp$X[, rand_var_names])
    model_data = list(N = nrow(tmp$X),
                      n_vars = ncol(tmp$X) - length(rand_var_names),
                      n_rand_var = length(rand_var_names),
                      n_rand_total = max(unlist(rand_idx)),
                      rand_var_shared = get_shared_rand(tmp$X[,rand_var_names]),
                      X = tmp$X[, !names(tmp$X) %in% rand_var_names],
                      rand_id = rand_idx,
                      Cq = tmp$Cq)
    return(model_data)
}

get_beta_intervals = function(fit, probs = c(0.025, 0.975),
                              pars = "betas")
{
    pp = summary(fit, probs = probs, pars = pars)$summary
    pp[,grep(".*%$", colnames(pp))]
}

est_power = function(mod, n_experiments, vars_list, betas, sigma_rand, sigma_Cq,
                     mod_matrix = expand.grid(vars_list),
                     rand_var_names = c("tech_rep", "rep"),
                     result_type = "power", # power or interval
                     iter = 400, probs = c(0.025, 0.975), ...)
{
    ans = replicate(n_experiments, {
        simulated_data = prep_model_data(, betas,
                                         sigma_rand, sigma_Cq,
                                         mod_matrix, rand_var_names)
        fit = sampling(mod, data = simulated_data, iter = iter)
        int = get_beta_intervals(fit, probs = probs)
        switch(result_type,
               power = apply(int, 1, excludes_zero),
               interval = sapply(seq_along(betas), function(i)
                   in_interval(int[i,], betas[i])))
    })
    power_calc(ans)
}

