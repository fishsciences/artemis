# Helper function to arrange data for Stan

prep_data = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas,
                     Cq_upper = 40, rand_sd = double(0),
                     prior_int, prior_b,
                     type = c("model", "sim"))
{
    has_inter = has_intercept(mod_list$x)
    x = remove_intercept(mod_list$x)
    n_vars = if(is.null(ncol(x))) 0 else ncol(x)
    priors = prep_priors(prior_b, x, mod_list$y)
    model_data = list(N = length(mod_list$y),
                      n_vars = n_vars,
                      X = x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper,
                      rand_sigma = as.array(rand_sd),
                      prior_mu = priors$location,
                      prior_sd = priors$scale,
                      has_inter = has_inter)

    model_data$prior_int_mu = prior_int$location
    model_data$prior_int_sd = prior_int$scale
    model_data$has_prior = 1L # for testing ifelse(length(b_prior_mu), 1, 0)
    
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

remove_intercept = function(x)
{
    cnms = colnames(x)
    i = grep("(Intercept)", cnms, invert = TRUE)

    x = as.data.frame(x[,i])
    colnames(x) = cnms[i]
    return(x)
}

prep_priors = function(prior, x, y)
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
