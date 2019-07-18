## These functions are for generating simulated data to be used
## to calculate power for the eDNA experiments
## M. Espe
## March 2019

sim_eDNA_lm = function(formula, vars_list,
                       betas, sigma_Cq,
                       std_curve_alpha, std_curve_beta,
                       n_sim = 1L,
                       upper_Cq = 40,
                       X = expand.grid(vars_list))
{
    if(!has_response(formula))
        stop("Please provide a dummy response variable for simulations")
    ml = gen_model_list_lm(formula, X)
    
    if(ncol(ml$x) != length(betas))
        stop("Please provide one beta per model term \n",
             "Provided: ", length(betas), "\n",
             "Required: ", ncol(ml$x), "\n")

    md = prep_sim(ml, std_curve_alpha, std_curve_beta, sigma_Cq, betas)
    sims = sampling(stanmodels$eDNA_sim, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim)

    sims = extract(sims)
    
    return(sims)
}

sim_eDNA_lmer = function(formula, vars_list,
                         betas, sigma_Cq,
                         sigma_rand,
                         std_curve_alpha, std_curve_beta,
                         n_sim = 1L,
                         upper_Cq = 40,
                         X = expand.grid(vars_list))
{
    if(!has_response(formula))
        stop("Please provide a dummy response variable for simulations")
    
    if(is_lme4(formula) && length(sigma_rand) == 0)
        stop("You must provide the sd of the random effects")

    ml = gen_model_list_lmer(formula, X)

    if(ncol(ml$x) != length(betas))
        stop("Please provide one beta per model term \n",
             "Provided: ", length(betas), "\n",
             "Required: ", ncol(ml$x), "\n")

    if(length(sigma_rand) != ncol(ml$groups))
        stop("You must provide one sd for each random effect\n",
             "Provided: ", length(sigma_rand), "\n",
             "Required: ", ncol(ml$groups), "\n",
             "Random effects: ", colnames(ml$groups))

    
    md = prep_sim(ml, std_curve_alpha, std_curve_beta, sigma_Cq,
                  betas = betas, rand_sd = sigma_rand)

    sims = sampling(stanmodels$eDNA_sim, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim, warmup = 0L)

    sims = extract(sims)
    
    return(sims)
}


gen_rand_eff = function(X, sigma_rand_eff)
{
    sapply(seq_along(X), function(i){
        n_levels = length(unique(X[[i]]))
        rnorm(n_levels, mean = 0, sd = sigma_rand_eff)[X[[i]]]
    })
}


gen_eDNA_conc = function(X, betas, rands)
    ## Generate the concentration from a
    ## model matrix X and a vector betas
{
    y_hat = as.matrix(X) %*% betas + rowSums(rands)
    return(y_hat)
}

gen_Cq = function(ln_conc, sigma,
                  std_alpha, std_beta,
                  thresh)
    ## Generate Cq from a conc, and then truncate to thresh
{
    Cq_hat = ln_std_curve(exp(ln_conc), std_alpha, std_beta)
    Cq_star = rnorm(length(ln_conc), Cq_hat, sigma)
    Cq_star[Cq_star > thresh] = thresh
    return(Cq_star)
}

ln_std_curve = function(x, alpha = 21.167769, beta = -1.52868305)
    ## equation, alpha, and beta from 20180927_VS_standard_curve.xlsx
{
    beta * log(x) + alpha
}

has_response = function(formula)
{
    attr(terms(formula), "response") > 0
}

prep_sim = function(mod_list, alpha, beta, Cq_sd, betas, rand_sd = double(0))
{
    model_data = list(N = length(mod_list$y),
                      n_vars = ncol(mod_list$x),
                      X = mod_list$x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      sigma_Cq = Cq_sd,
                      betas = betas)

    if(is.null(mod_list$groups)){
        b = list(model_data, 
                 has_rand = 0L,
                 n_rand_var = 0L,
                 n_rand_total = 0L, 
                 rand_var_shared = double(0),
                 groups = double(0),
                 rand_sigma = rand_sd)
        model_data = c(model_data, b)
    } else { 
        rand_idx = unlist(relevel_rands(mod_list$groups))
        rand_var_shared = get_shared_rand(mod_list$groups)
        
        b = list(has_rand = 1L,
                 n_rand_var = ncol(mod_list$groups),
                 n_rand_total = max(rand_idx),
                 rand_var_shared = rand_var_shared,
                 groups = rand_idx,
                 rand_sigma = as.array(rand_sd))
        
        model_data = c(model_data, b)
    }

    return(model_data)
}
