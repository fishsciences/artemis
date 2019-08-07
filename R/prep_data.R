
prep_data = function(mod_list,
                     alpha, beta,
                     Cq_sd, betas,
                     Cq_upper = 40, rand_sd = double(0),
                     type = c("model", "sim"))
{
    model_data = list(N = length(mod_list$y),
                      n_vars = ncol(mod_list$x),
                      X = mod_list$x,
                      std_curve_alpha = alpha,
                      std_curve_beta = beta,
                      upper_Cq = Cq_upper,
                      rand_sigma = as.array(rand_sd))

    if(is.null(mod_list$groups)){
        model_data$has_rand = 0L
        model_data$n_rand_var = 0L
        model_data$n_rand_total = 0L
        model_data$rand_var_shared = double(0)
        model_data$groups = double(0)
    } else { 
        rand_idx = unlist(relevel_rands(mod_list$groups))
        rand_var_shared = get_shared_rand(mod_list$groups)
        
        model_data$n_rand = mod_list$n_rand
        model_data$n_rand_total = max(rand_idx)
        model_data$rand_var_shared = rand_var_shared
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
