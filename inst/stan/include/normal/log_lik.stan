  for (n in 1:N_obs)
    log_lik[n] = normal_lpdf(y_obs[n] | mu_all[n], sigma_ln_eDNA);
  for (n in 1:N_cens)
    log_lik[N_obs + n] = normal_lcdf(L[n] | mu_all[N_obs + n], sigma_ln_eDNA);
