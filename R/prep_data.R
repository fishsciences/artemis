# Helper function to arrange data for Stan

prep_data = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas,
                     Cq_upper = 40, rand_sd = double(0),
                     b_prior_mu, b_prior_sd,
                     type = c("model", "sim"))
{
    has_inter = has_intercept(mod_list$x)
    x = remove_intercept(mod_list$x)
    n_vars = if(is.null(ncol(x))) 0 else ncol(x)
    model_data = list(N = length(mod_list$y),
                      n_vars = n_vars,
                      X = x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper,
                      rand_sigma = as.array(rand_sd),
                      prior_mu = b_prior_mu,
                      prior_sd = b_prior_sd,
                      has_inter = has_inter)

    ## Temp for testing
    model_data$prior_int_mu = 0
    model_data$prior_int_sd = 10
    model_data$has_prior = ifelse(length(b_prior_mu), 1, 0)
    model_data$n_prior = length(b_prior_mu)
    
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
