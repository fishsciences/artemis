
est_p_detect = function(var_levels,
                        betas, 
                        Cq_sd,
                        std_curve_alpha, std_curve_beta,
                        n_rep = 1:12,
                        model_fit = NULL,
                        target_detection = 0.8,
                        upper_Cq = 40)

{
    if(!is.null(dim(var_levels)))
        stop("Sorry, only one set of variable levels at a time currently supported")
    if(is.null(model_fit) && length(var_levels) != length(betas))
        stop("Variable levels and betas cannot be of different lengths")

    ln_conc_hat = if(is.null(model_fit)) {
                      var_levels %*% betas
                  } else {
                      apply(model_fit@betas, 1, function(y) var_levels %*% y)
                  }

    ln_thresh = (upper_Cq - std_curve_alpha) / std_curve_beta

    ans = sapply(n_rep, function(i) 1 - (pnorm(ln_thresh, ln_conc_hat, Cq_sd) ^ i))

    structure(ans,
              variable_levels = var_levels,
              reps = n_rep,
              class = c("eDNA_p_detect", class(ans)))
}
