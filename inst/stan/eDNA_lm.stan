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

transformed data {
  matrix[N_obs, K] Q_ast_obs;
  matrix[K, K] R_ast_obs;
  matrix[K, K] R_ast_inverse_obs;
  matrix[N_cens, K] Q_ast_cens;
  matrix[K, K] R_ast_cens;
  matrix[K, K] R_ast_inverse_cens;
  // thin and scale the QR decomposition
  Q_ast_obs = qr_thin_Q(X_obs) * sqrt(N_obs - 1);
  R_ast_obs = qr_thin_R(X_obs) / sqrt(N_obs - 1);
  R_ast_inverse_obs = inverse(R_ast_obs);
  Q_ast_cens = qr_thin_Q(X_cens) * sqrt(N_cens - 1);
  R_ast_cens = qr_thin_R(X_cens) / sqrt(N_cens - 1);
  R_ast_inverse_cens = inverse(R_ast_cens);
}
parameters {
  real intercept[has_inter ? 1 : 0];
  vector[K] thetas;      // coefficients on Q_ast
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
  vector[K] betas;
  betas = R_ast_inverse_obs * thetas; // coefficients on x
}

model {
  vector[N_obs] mu_obs = (K ? Q_ast_obs * thetas : rep_vector(0.0, N_obs)) +
	(has_inter ? intercept[1] : 0.0);
  vector[N_cens] mu_cens = (K ? Q_ast_cens * thetas : rep_vector(0.0, N_cens)) +
	(has_inter ? intercept[1] : 0.0);

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);

  for(k in 1:K)
	betas[k] ~ normal(prior_mu[k], prior_sd[k]);

  sigma_ln_eDNA ~ normal(0, 1);
  
  target += normal_lpdf(y_obs | mu_obs, sigma_ln_eDNA);
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
