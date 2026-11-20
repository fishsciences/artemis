/* Stan model for a Poisson distributed linear mixed model
   M. Espe
*/

functions {
#include functions/zero_sum_components.stanfunctions
}

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
#include common/prior_data.stan
#include random/data.stan
  matrix[N, K_r] X_r;
  array[N] int<lower=0> y;
}

parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] betas;
#include random/params.stan
}

transformed parameters {
#include random/tparams.stan
}

model {
  // priors
#include common/fixed_priors.stan
#include random/priors.stan

  // Likelihood
  y ~ poisson_log_glm(append_col(X, X_r),
                      (has_inter ? intercept[1] : 0.0),
                      append_row(betas, rand_betas));
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] |
                                  append_col(X[n], X_r[n]) *
                                  append_row(betas, rand_betas) +
                                  (has_inter ? intercept[1] : 0.0));
}
