/*
  Censored linear model on ln[eDNA]
 */

data {
#include normal/data.stan
#include common/prior_data.stan
}

transformed data {
#include normal/tdata_qr.stan
}

parameters {
#include normal/params.stan
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
#include normal/tparams.stan
}

model {
#include normal/mu.stan

  // priors
#include common/fixed_priors.stan
  sigma_ln_eDNA ~ normal(0, 1);

#include normal/lik.stan
}

generated quantities {
  vector[N_obs + N_cens] log_lik;
  {
#include normal/mu.stan
#include normal/log_lik.stan
  }
}
