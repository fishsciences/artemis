# These functions are for generating simulated data to be used
# to calculate power for the eDNA experiments
# M. Espe
# March 2019

sim_eDNA_lm = function(formula, vars_list,
                       betas, sigma_Cq,
                       std_curve_alpha, std_curve_beta,
                       upper_Cq = 40,
                       X = expand.grid(vars_list))
{
    ml = gen_model_list_lm(formula, X)

    if(ncol(ml$x) != length(betas))
        stop("Please provide one beta per model term")
    
    ln_conc_hat = ml$x %*% betas

    cq_star = gen_Cq(ln_conc_hat, sigma_Cq, 
                     std_curve_alpha, std_curve_beta,
                     thresh = upper_Cq)
    
    return(list(x = X, Cq_star = cq_star, ln_conc = ln_conc_hat))

}

sim_eDNA_lmer = function(formula, vars_list,
                       X = expand.grid(vars_list),
                       betas, sigma_Cq,
                       sigma_rand = numeric())
{
    if(is_lme4(formula) && length(sigma_rand) == 0)
        stop("You must provide the sd of the random effects")
}

gen_rand_eff = function(X, sigma_rand_eff)
{
    sapply(seq_along(X), function(i){
        n_levels = length(unique(X[[i]]))
        rnorm(n_levels, mean = 0, sd = sigma_rand_eff)[X[[i]]]
    })
}


gen_eDNA_conc = function(X, betas, rands)
    # Generate the concentration from a
    # model matrix X and a vector betas
{
    y_hat = as.matrix(X) %*% betas + rowSums(rands)
    return(y_hat)
}

gen_Cq = function(ln_conc, sigma,
                  std_alpha, std_beta,
                  thresh)
    # Generate Cq from a conc, and then truncate to thresh
{
    Cq_hat = ln_std_curve(exp(ln_conc), std_alpha, std_beta)
    Cq_star = rnorm(length(ln_conc), Cq_hat, sigma)
    Cq_star[Cq_star > thresh] = thresh
    return(Cq_star)
}

ln_std_curve = function(x, alpha = 21.167769, beta = -1.52868305)
    # equation, alpha, and beta from 20180927_VS_standard_curve.xlsx
{
    beta * log(x) + alpha
}
    
