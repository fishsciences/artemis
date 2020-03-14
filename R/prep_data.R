# Helper function to arrange data for Stan

prep_data = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas,
                     Cq_upper = 40, rand_sd = double(0),
                     b_prior_mu, b_prior_sd,
                     type = c("model", "sim"))
{
    model_data = list(N = length(mod_list$y),
                      n_vars = ncol(mod_list$x),
                      X = mod_list$x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper,
                      rand_sigma = as.array(rand_sd),
                      prior_mu = b_prior_mu,
                      prior_sd = b_prior_sd)

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
