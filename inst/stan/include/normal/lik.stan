  target += normal_lpdf(y_obs | mu_all[1:N_obs], sigma_ln_eDNA);
  if (N_cens)
    target += normal_lcdf(L | mu_all[(N_obs + 1):], sigma_ln_eDNA);
