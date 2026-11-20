  {
    real nz_a = has_nz_inter ? nz_alpha[1] : 0.0;
    target += bernoulli_logit_glm_lpmf(1 | X_nz, nz_a, nz_beta) +
      normal_lpdf(y_obs | mu_all[1:N_obs], sigma_ln_eDNA);
    if (N_cens)
      target += log_sum_exp(bernoulli_logit_glm_lpmf(0 | X_z, nz_a, nz_beta),
                            bernoulli_logit_glm_lpmf(1 | X_z, nz_a, nz_beta) +
                            normal_lcdf(L | mu_all[(N_obs + 1):], sigma_ln_eDNA));
  }
