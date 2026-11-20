/*
  Censored linear mixed model on ln[eDNA]
 */

functions {
#include functions/num_nonzero.stanfunctions
#include functions/zero_sum_components.stanfunctions
}

data {
#include normal/data.stan
#include random/data.stan
  matrix[N_obs, K_r] X_obs_r;
  matrix[N_cens, K_r] X_cens_r;
#include common/prior_data.stan
}

transformed data {
#include random/tdata_csr.stan
#include normal/tdata_qr.stan
}

parameters {
#include normal/params.stan
#include random/params.stan
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
#include normal/tparams.stan
#include random/tparams.stan
}

model {
#include normal/mu.stan
#include random/mu.stan

  // priors
#include common/fixed_priors.stan
#include random/priors.stan
  sigma_ln_eDNA ~ exponential(1);

#include normal/lik.stan
}

generated quantities {
  vector[N_obs + N_cens] log_lik;
  {
#include normal/mu.stan
#include random/mu.stan
#include normal/log_lik.stan
  }
}
