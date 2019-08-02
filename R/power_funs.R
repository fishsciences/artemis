
est_p_detect = function(var_levels,
                        betas, 
                        Cq_sd,
                        std_curve_alpha, std_curve_beta,
                        n_rep = 1:12,
                        model_fit = NULL,
                        target_detection = 0.8,
                        upper_Cq = 40,
                        n_sim = 500L)

{
    ln_conc_hat = if(is.null(model_fit)) {
                      var_levels %*% betas
                  } else {
                      apply(model_fit@betas, 1, function(y) var_levels %*% y)
                  }

    ln_thresh = (upper_Cq - std_curve_alpha) / std_curve_beta


    sapply(n_rep, function(i) 1 - (pnorm(ln_thresh, ln_conc_hat, Cq_sd) ^ i))
}
