  intercept ~ normal(prior_int_mu, prior_int_sd);
  if (K) {
    for (k in 1:K)
      betas[k] ~ normal(prior_mu[k], prior_sd[k]);
  }
