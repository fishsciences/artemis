/* Stan model for a Poisson distributed linear model
   M. Espe
*/

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, K] X;
#include common/prior_data.stan
  array[N] int<lower=0> y;
}

parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] betas;
}

model {
  // priors
#include common/fixed_priors.stan

  // Likelihood
  y ~ poisson_log_glm(X, (has_inter ? intercept[1] : 0.0), betas);
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N)
    log_lik[n] = poisson_log_lpmf(y[n] | X[n] * betas + (has_inter ? intercept[1] : 0.0));
}
