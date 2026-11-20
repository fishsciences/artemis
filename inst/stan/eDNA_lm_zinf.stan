/*
  Zero-inflated censored linear model on ln[eDNA]
 */

data {
#include normal/data.stan
#include zinf/data.stan
#include common/prior_data.stan
}

transformed data {
#include normal/tdata_qr.stan
}

parameters {
#include normal/params.stan
#include zinf/params.stan
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
#include normal/tparams.stan
}

model {
#include normal/mu.stan

  // priors
#include common/fixed_priors.stan
#include zinf/priors.stan
  sigma_ln_eDNA ~ std_normal();

#include zinf/lik.stan
}

generated quantities {
  vector[N_obs + N_cens] log_lik;
  {
#include normal/mu.stan
#include zinf/log_lik.stan
  }
}
