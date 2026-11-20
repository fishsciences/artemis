  {
    real nz_a = has_nz_inter ? nz_alpha[1] : 0.0;
    for (n in 1:N_obs)
      log_lik[n] = bernoulli_logit_glm_lpmf(1 | to_matrix(X_nz[n]), nz_a, nz_beta) +
        normal_lpdf(y_obs[n] | mu_all[n], sigma_ln_eDNA);
    for (n in 1:N_cens)
      log_lik[N_obs + n] =
        log_sum_exp(bernoulli_logit_glm_lpmf(0 | to_matrix(X_z[n]), nz_a, nz_beta),
                    bernoulli_logit_glm_lpmf(1 | to_matrix(X_z[n]), nz_a, nz_beta) +
                    normal_lcdf(L[n] | mu_all[N_obs + n], sigma_ln_eDNA));
  }
