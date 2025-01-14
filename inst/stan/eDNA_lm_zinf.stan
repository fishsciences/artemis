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

  // Zero-inf components
  int<lower=0> K_z; // predictors
  matrix[N_cens, K_z] X_z;
  matrix[N_obs, K_z] X_nz;
  
  // prior parameters - user provided
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;  

  // for intercept-less models
  int<lower=0,upper=1> has_inter;
}

transformed data {
  matrix[N_obs + N_cens, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  // thin and scale the QR decomposition
  if(K){
	Q_ast = qr_thin_Q(append_row(X_obs, X_cens)) * sqrt((N_obs+N_cens) - 1);
	R_ast = qr_thin_R(append_row(X_obs, X_cens)) / sqrt((N_obs+N_cens) - 1);
	R_ast_inverse = inverse(R_ast);
  }
}
parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] thetas;      // coefficients on Q_ast
  real<lower=0, upper=1> z_alpha;
  vector[K_z] z_beta;
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
  vector[K] betas;
  if(K){
    betas = R_ast_inverse * thetas; // coefficients on x
  }
}

model {
  vector[N_obs + N_cens] mu_all = (K ? Q_ast * thetas : rep_vector(0.0, N_obs + N_cens)) +
    (has_inter ? intercept[1] : 0.0);
  vector[N_obs] mu_obs = mu_all[1:N_obs];
  vector[N_cens ? N_cens : 1] mu_cens = N_cens ? mu_all[(N_obs+1):] : [0]';

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);

  if(K) {
    for(k in 1:K)
      betas[k] ~ normal(prior_mu[k], prior_sd[k]);
  }
  
  sigma_ln_eDNA ~ normal(0, 1);

  // non-zero - efficient vectorized
  target += bernoulli_logit_glm_lpmf(1 | X_nz, z_alpha, z_beta) +
	normal_lpdf(y_obs | mu_obs, sigma_ln_eDNA);

  // zero - inefficient looping
  if(N_cens){
	array[1] int zeros = rep_array(0,1);
	array[1] int ones = rep_array(1,1);
	for (n in 1:N_cens){
	  target += log_sum_exp(bernoulli_logit_glm_lpmf(zeros[1] | to_matrix(X_z[n]), z_alpha, z_beta),
							bernoulli_logit_glm_lpmf(ones[1] | to_matrix(X_z[n]), z_alpha, z_beta) +
							normal_lcdf(L | mu_cens[n], sigma_ln_eDNA));
	}
  }
}

generated quantities{
  vector[N_obs + N_cens] log_lik = rep_vector(0, N_obs + N_cens);

  if(N_obs){
	array[1] int zeros = rep_array(0,1);
	array[1] int ones = rep_array(1,1);
    for(n in 1:N_obs){
      log_lik[n] = bernoulli_logit_glm_lpmf(ones[1] | to_matrix(X_nz[n]), z_alpha, z_beta) +
		normal_lpdf(y_obs[n] | (has_inter ? intercept[1] : 0.0) +
					(K ? X_obs[n] * betas : 0), sigma_ln_eDNA);
    }
  }
  if(N_cens){
   	array[1] int zeros = rep_array(0,1);
	array[1] int ones = rep_array(1,1);
	for(n in 1:N_cens){
      log_lik[n+N_obs] = log_sum_exp(bernoulli_logit_glm_lpmf(zeros[1] | to_matrix(X_z[n]), z_alpha, z_beta),
									 bernoulli_logit_glm_lpmf(ones[1] | to_matrix(X_z[n]), z_alpha, z_beta) +
									 normal_lcdf(L[n] | (has_inter ? intercept[1] : 0) +
												 (K ? X_cens[n] * betas : 0), sigma_ln_eDNA));
    }
  }
}
