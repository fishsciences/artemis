/*
  New, simplified model
 */
functions{
  int num_nonzero(matrix m){
	int n_row = rows(m);
	int n_col = cols(m);
	int ans = 0;

	for(i in 1:n_row)
	  for(j in 1:n_col)
		if(m[i,j] > 0)
		  ans += 1;
	
	return ans;
  }

}

data {
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  int<lower=0> K; // predictors
  matrix[N_obs, K] X_obs;
  matrix[N_cens, K] X_cens;
  vector[N_obs] y_obs; // ln[eDNA]
  vector[N_cens] L; // lower bound on ln[eDNA]
  //#include lm_data.stan // Play with this later
  // Random effects
  int<lower=0> K_r;
  int<lower=1> N_grp;
  int<lower=1,upper=N_grp> group[K_r];
  matrix[N_obs, K_r] X_obs_r;
  matrix[N_cens, K_r] X_cens_r;
  
  // prior parameters - user provided
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;  
  real rand_sd;
  
  // for intercept-less models
  int<lower=0,upper=1> has_inter;
}

transformed data {
  // QR decomp
  matrix[N_obs + N_cens, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;

  // Random effect X to sparse mat for efficiency
  matrix[N_obs+N_cens, K_r] X_r_all = append_row(X_obs_r, X_cens_r);
  int n_nonzero = num_nonzero(X_r_all);
  vector[n_nonzero] w = csr_extract_w(X_r_all);
  int v[n_nonzero] = csr_extract_v(X_r_all);
  int u[N_obs+N_cens+1] = csr_extract_u(X_r_all);

  if(K){
	// thin and scale the QR decomposition
	Q_ast = qr_thin_Q(append_row(X_obs, X_cens)) * sqrt((N_obs+N_cens) - 1);
	R_ast = qr_thin_R(append_row(X_obs, X_cens)) / sqrt((N_obs+N_cens) - 1);
	R_ast_inverse = inverse(R_ast);
  }
}
parameters {
  real intercept[has_inter ? 1 : 0];
  vector[K] thetas; // coefficients on Q_ast
  vector[K_r] rand_betas_raw;
  vector<lower=0>[N_grp] rand_sigma;
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
  vector[K] betas;
  vector[K_r] rand_betas;
  if(K)
	betas = R_ast_inverse * thetas; // coefficients on x
  for(k in 1:K_r)
	rand_betas[k] = rand_betas_raw[k] * rand_sigma[group[k]];
}

model {
  vector[N_obs + N_cens] mu_all =
	(K ? Q_ast * thetas : rep_vector(0.0, N_obs + N_cens)) +
	csr_matrix_times_vector(N_obs+N_cens, K_r, w, v, u, rand_betas) +
  /* (append_row(X_obs_r, X_cens_r) * rand_betas) + */
	(has_inter ? intercept[1] : 0.0);
  vector[N_obs] mu_obs = mu_all[1:N_obs];
  vector[N_cens] mu_cens = mu_all[(N_obs+1):];

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);
  rand_sigma ~ gamma(2, 0.1);

  if(K)
	for(k in 1:K)
	  betas[k] ~ normal(prior_mu[k], prior_sd[k]);

  rand_betas_raw ~ std_normal();
  
  sigma_ln_eDNA ~ exponential(1);
  
  target += normal_lpdf(y_obs | mu_obs, sigma_ln_eDNA);
  target += normal_lcdf(L | mu_cens, sigma_ln_eDNA);
}

generated quantities{
  vector[N_obs + N_cens] log_lik;
  {
    matrix[N_obs, K + K_r] X_obs_tmp = append_col(X_obs, X_obs_r);
    matrix[N_cens, K + K_r] X_cens_tmp = append_col(X_cens, X_cens_r);
    vector[K+K_r] betas_tmp = append_row(betas, rand_betas);
    
    for(n in 1:N_obs){
      
      log_lik[n] = normal_lpdf(y_obs[n] | (has_inter ? intercept[1] : 0.0) +
			       X_obs_tmp[n] * betas_tmp, sigma_ln_eDNA);
    }
    for(n in 1:N_cens){
      log_lik[n+N_obs] = normal_lcdf(L[n] | (has_inter ? intercept[1] : 0) +
				     X_cens_tmp[n] * betas_tmp, sigma_ln_eDNA);
    }
  }
}
