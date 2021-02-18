/*
  New, simplified model
 
  Zero-inflated
*/

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
  vector[K] prior_mu;
  vector[K] prior_sd;  

  
}
parameters {
  real intercept;
  vector[K] betas;
  real<lower=0> sigma_Cq;
  real<lower=0, upper=1> p_zero;
}
model {
  /* vector[N_obs] mu_obs = intercept + X_obs * betas; */
  vector[N_cens] mu_cens = intercept + X_cens * betas;

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);

  for(k in 1:K)
	betas[k] ~ normal(prior_mu[k], prior_sd[k]);

  sigma_Cq ~ normal(0, 1);
  
  y_obs ~ normal_id_glm(X_obs, intercept, betas, sigma_Cq);
  // zero-inflated
  target += N_obs * bernoulli_lpmf(0 | p_zero);
  for(n in 1:N_cens)
	target += log_sum_exp(bernoulli_lpmf(1 | p_zero),
						  bernoulli_lpmf(0 | p_zero) +
						  normal_lcdf(L | mu_cens[n], sigma_Cq));
}

generated quantities{
  vector[N_obs + N_cens] log_lik;

  for(n in 1:N_obs){
	log_lik[n] = normal_lpdf(y_obs[n] | intercept + X_obs[n] * betas, sigma_Cq);
  }
  
  for(n in 1:N_cens){
	log_lik[n+N_obs] = log_sum_exp(bernoulli_lpmf(1 | p_zero),
								   bernoulli_lpmf(0 | p_zero) +
								   normal_lcdf(L | intercept + X_cens[n] * betas,
											   sigma_Cq));
	
  }
}

