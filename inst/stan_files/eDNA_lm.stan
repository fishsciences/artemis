data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K; // predictors
  matrix[N_obs, K] X_obs;
  matrix[N_cens, K] X_cens;
  vector[N_obs] y_obs; // ln[eDNA]
  real<upper=min(y_obs)> L; // lower bound on ln[eDNA]

  // prior parameters - user provided
  real prior_int_mu;
  real prior_int_sd;
}
parameters {
  real intercept;
  vector[K] betas;
  real<lower=0> sigma_Cq;
}
model {
  /* vector[N_obs] mu_obs = intercept + X_obs * betas; */
  vector[N_cens] mu_cens = intercept + X_cens * betas;

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);
  betas ~ normal(0, 10);
  sigma_Cq ~ std_normal();
  
  y_obs ~ normal_id_glm(X_obs, intercept, betas, sigma_Cq);
  target += normal_lcdf(L | mu_cens, sigma_Cq);
}
