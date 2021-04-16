## Helper function to arrange data for Stan
#The mode and sim functions
# have diverged to the point that the copied code is worth having the
# sim and model prep functions separate from eachother and cleaner.

prep_data.lmer = function(mod_list,
                          alpha, beta,
                          prob_zero = 0,
                          Cq_upper,
                          prior_int, prior_b,
                          rand_sd)
{
    # This gets most of the data in there
    model_data = prep_data.lm(mod_list, alpha, beta, 
                             prob_zero, Cq_upper, prior_int, prior_b)
    i = mod_list$y < Cq_upper

    model_data$X_obs_r = mod_list$rand_x[i,]
    model_data$X_cens_r = mod_list$rand_x[!i,]

    model_data$K_r = mod_list$n_rand
    model_data$group = mod_list$groups
    model_data$N_grp = mod_list$n_grp
    if(rand_sd$dist != "exponential")
        stop("Only exponential priors on random effect variance supported at this time")
    model_data$rand_sd = rand_sd$scale^-1    
    return(model_data)
}

# This needs to be updated maybe?
prep_data.sim = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas, p_zero,
                     Cq_upper = 40, rand_sd = double(0),
                     prior_int, prior_b, error_type = "fixed")
{
        
    has_inter = has_intercept(mod_list$x)
    x = mod_list$x
    idx = seq(length(mod_list$y))
    n_below = 1
    
    n_vars = if(is.null(ncol(x))) 0 else ncol(x)
    priors = prep_priors(prior_b, x, mod_list$y)
    
    model_data = list(y = mod_list$y,
                      N = length(mod_list$y),
                      n_vars = n_vars,
                      X = as.matrix(x[idx,]),
                      upper_Cq = Cq_upper,
                      p_zero = p_zero,
                      rand_sigma = as.array(rand_sd),
                      prior_mu = priors$location,
                      prior_sd = priors$scale,
                      has_inter = has_inter,
                      n_below = n_below)

    if(length(alpha) < model_data$N)
        alpha = rep(alpha, model_data$N)
    if(length(beta) < model_data$N)
        beta = rep(beta, model_data$N)

    model_data$std_curve_alpha = alpha
    model_data$std_curve_beta = beta
    
    model_data$prior_int_mu = prior_int$location
    model_data$prior_int_sd = prior_int$scale
    model_data$has_prior = 1L # for testing ifelse(length(b_prior_mu), 1, 0)

    model_data$sd_vary = switch(error_type,
                                "fixed" = 0L,
                                "varying" = 1L,
                                stop("CQ Error type must be either 'fixed' or 'varying'"))
    model_data$center_Cq = 30 ## for testing
    
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
        model_data$rand_x = mod_list$rand_x[idx,]
        model_data$n_grp = mod_list$n_grp
        model_data$rand_sigma = as.array(rand_sd)
    }
        model_data$betas = as.array(betas)
        model_data$sigma_ln_eDNA = Cq_sd
    
    return(model_data)
}

remove_intercept = function(x)
{
    cnms = colnames(x)
    i = grep("(Intercept)", cnms, invert = TRUE)

    x = as.data.frame(subset(x, select = i)) ## Fixes dropping names on 1 column df
    
    return(x)
}

prep_priors = function(prior, x, y)
## Follows general advice from rstanarm package
{
    n = ncol(x)
    if(n == 0){
        prior$location = double(0)
        prior$scale = double(0)
        return(prior)
    }
    
    if(length(prior$location) < n)
        prior$location = rep(prior$location, n)
    if(length(prior$scale) < n)
        prior$scale = rep(prior$scale, n)

    if(prior$autoscale){
        prior$scale = prior$scale * sd(y)
        multi_val = apply(x, 2, function(x) length(unique) > 2)
        x_sd = apply(x, 2, sd)
        prior$scale[multi_val] = prior$scale[multi_val] * x_sd[multi_val] 
    }
    prior[c("location","scale")] = lapply(prior[c("location","scale")], as.array)
    
    return(prior)
}


prep_data.lm = function(mod_list,
                        alpha, beta,
                        prob_zero = 0,
                        Cq_upper,
                        prior_int, prior_b,
                        rand_sd = 0) # ignored for lm

{
    has_inter = has_intercept(mod_list$x)
    x = remove_intercept(mod_list$x)
    i = mod_list$y < Cq_upper
    
    n_vars = if(is.null(ncol(x))) 0 else ncol(x)
    
    # Convert to ln[eDNA] here to avoid thinking too hard in the model
    ln_eDNA = (mod_list$y - alpha) / beta
    priors = prep_priors(prior_b, x, ln_eDNA)

    lower_bound = (Cq_upper - alpha) / beta

    # to deal with multiple std curves, we need expand the lower bound
    # to length N
    if(length(lower_bound) == 1) lower_bound = rep(lower_bound, length(ln_eDNA))
    
    model_data = list(y_obs = ln_eDNA[i],
                      N_obs = sum(i),
                      N_cens = sum(!i),
                      K = n_vars,
                      X_obs = as.matrix(x[i, , drop = FALSE],),
                      X_cens = as.matrix(x[!i, , drop = FALSE],),
                      L = lower_bound[!i],
                      p_zero = prob_zero,
                      prior_mu = priors$location,
                      prior_sd = priors$scale,
                      has_inter = has_inter,
                      X = x)
    
    ## avoid #837 in rstan
    if(typeof(model_data$X_obs) == "logical" && ncol(model_data$X_obs) == 0)
        storage.mode(model_data$X_obs) = "numeric"
    if(typeof(model_data$X_cens) == "logical" && ncol(model_data$X_cens) == 0)
        storage.mode(model_data$X_cens) = "numeric"
    
    model_data$prior_int_mu = prior_int$location
    model_data$prior_int_sd = prior_int$scale

    model_data$center_Cq = 30 ## for testing
   
    return(model_data)
}
