# These functions are for generating simulated data to be used
# to calculate power for the eDNA experiments
# M. Espe
# March 2019

sim_data = function(formula, vars_list, betas,
                    sigma_rand, sigma_Cq,
                    X = expand.grid(vars_list),
                    rand_var_names = c("tech_rep", "rep"))
{
    if(ncol(X[,!names(X) %in% rand_var_names]) != length(betas))
        stop("'betas' must the same length as number of fixed effects")
    
    rand = gen_rand_eff(X[,rand_var_names], sigma_rand)
    ln_conc = gen_eDNA_conc(X[, !names(X) %in% rand_var_names],
                            betas, rand)
    Cq = gen_Cq(ln_conc, sigma_Cq)
    return(list(X = X, ln_conc = ln_conc, Cq = Cq))
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
                  std_curve_fun = ln_std_curve,
                  thresh = 40)
    # Generate Cq from a conc, and then truncate to thresh
{
    Cq_hat = std_curve_fun(exp(ln_conc))
    Cq_star = rnorm(length(ln_conc), Cq_hat, sigma)
    Cq_star[Cq_star > thresh] = thresh
    return(Cq_star)
}

ln_std_curve = function(x, alpha = 21.167769, beta = -1.52868305)
    # equation, alpha, and beta from 20180927_VS_standard_curve.xlsx
{
    beta * log(x) + alpha
}
