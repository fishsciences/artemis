/*
  New, simplified model
 */

data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K; // predictors
  matrix[N_obs, K] X_obs;
  matrix[N_cens, K] X_cens;
  vector[N_obs] y_obs; // ln[eDNA]
  vector[N_cens] L; // lower bound on ln[eDNA]

  // prior parameters - user provided
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;  

  // for intercept-less models
  int<lower=0,upper=1> has_inter;
}

parameters {
  real intercept[has_inter ? 1 : 0];
  vector[K] betas;
  real<lower=0> sigma_ln_eDNA;
}

model {
  /* vector[N_obs] mu_obs = intercept + X_obs * betas; */
  vector[N_cens] mu_cens = (K ? X_cens * betas : rep_vector(0.0, N_cens)) +
	(has_inter ? intercept[1] : 0.0);

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);

  for(k in 1:K)
	betas[k] ~ normal(prior_mu[k], prior_sd[k]);

  sigma_ln_eDNA ~ normal(0, 1);
  
  y_obs ~ normal_id_glm(X_obs, has_inter ? intercept[1] : 0.0, betas, sigma_ln_eDNA);
  target += normal_lcdf(L | mu_cens, sigma_ln_eDNA);
}

generated quantities{
  vector[N_obs + N_cens] log_lik;

  for(n in 1:N_obs){
	log_lik[n] = normal_lpdf(y_obs[n] | (has_inter ? intercept[1] : 0.0) +
							 X_obs[n] * betas, sigma_ln_eDNA);
  }
  for(n in 1:N_cens){
	log_lik[n+N_obs] = normal_lcdf(L[n] | (has_inter ? intercept[1] : 0) +
								   X_cens[n] * betas, sigma_ln_eDNA);
  }
}
