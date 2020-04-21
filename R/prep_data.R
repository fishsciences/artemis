# Helper function to arrange data for Stan

prep_data = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas,
                     Cq_upper = 40, rand_sd = double(0),
                     prior_int, prior_beta, qr, 
                     type = c("model", "sim"))
{
    model_data = list(N = length(mod_list$y),
                      n_vars = ncol(mod_list$x),
                      X = mod_list$x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper,
                      rand_sigma = as.array(rand_sd),
                      intercept_mu = prior_int$location,
                      intercept_sd = prior_int$scale)
    prior_beta = prep_priors(prior_beta, mod_list$x)
    model_data$prior_mu = prior_beta$location
    model_data$prior_sd = prior_beta$scale
    model_data$has_intercept = has_intercept(mod_list$x)
    
    ## placeholder for HS priors
    model_data[c("global_scale", "nu_global","nu_local", "slab_scale", "slab_df")] = 1.0
    model_data$hs_prior = 0L
    
    model_data$has_prior = as.integer(!qr)
            
    if(is.null(mod_list$groups)){
        model_data$has_random = 0L
        model_data$n_rand = 0L
        model_data$groups = as.array(integer(0))
        model_data$rand_x = array(double(0), dim = c(nrow(mod_list$x), 0))
        model_data$n_grp = 0L
        model_data$rand_sigma = as.array(rand_sd)
    } else { 
        model_data$has_random = 1L
        model_data$n_rand = mod_list$n_rand
        model_data$groups = mod_list$groups
        model_data$rand_x = mod_list$rand_x
        model_data$n_grp = mod_list$n_grp
        model_data$rand_sigma = as.array(rand_sd)
    }
    if(type == "sim"){
        model_data$betas = as.array(betas)
    }

    return(model_data)
}


prep_priors = function(priors, X)
{
    n = ncol(X)
    
    if(!all(c("location", "scale", "df") %in% names(priors)))
        stop("Priors must be a named list with elements 'location', 'scale', and 'df'")

    if(length(priors$location) != 1 && length(priors$location) != n)
        stop("Prior location must either be length 1 or length ", n)

    if(length(priors$location) != n)
        priors$location = rep(priors$location, n)

    if(length(priors$scale) != n)
        priors$scale = rep(priors$scale, n)
    
    return(priors)
}
