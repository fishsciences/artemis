/**
 * New, simplified model
 */

// This needs to be copied into each file because there is a bug in the #includes parsing 

functions {
  int num_nonzero(matrix m){
	int n_row = rows(m);
	int n_col = cols(m);
	int ans = 0;
	
	for(i in 1:n_row)
	  for(j in 1:n_col)
		if(m[i , j] > 0)
		  ans += 1;
	
	return ans;
  }
    /**
   * Get the component size
   * Assumes the componenets are contiguous
   * required for segmented sum to zero
   */

  array[] int component_size(array[] int grp){
	int N = size(grp);
	int K = max(grp);
	array[K] int ans = zeros_int_array(K);
	
	for(n in 1:N){
	  ans[grp[n]] += 1;
	}

	return ans;
  }

  // Adapted from https://github.com/mitzimorris/geomed_2024/blob/main/stan/bym2_multicomp.stan
   /**
   * Component-wise constrain sum-to-zero vectors
   *
   * @param phi unconstrained vector of zero-sum slices
   * @param sizes component sizes
   * @return vector phi_ozs, the vector whose slices sum to zero
   */
  vector zero_sum_components(vector phi, array[] int sizes) {
    vector[sum(sizes)] phi_ozs;
    int idx_phi = 1;
    int idx_ozs = 1;
    for (i in 1:size(sizes)) {
      int n = sizes[i];
      phi_ozs[idx_ozs : (idx_ozs + n - 1)] = 
        zero_sum_constrain(segment(phi, idx_phi, n - 1));
      idx_phi += n - 1;
      idx_ozs += n;
    }
    return phi_ozs;
  }

  // From https://github.com/mitzimorris/geomed_2024/blob/main/stan/bym2_multicomp.stan
  /**
   * Constrain sum-to-zero vector
   *
   * @param y unconstrained zero-sum parameters
   * @return vector z, the vector whose slices sum to zero
   */
  vector zero_sum_constrain(vector y) {
    int N = num_elements(y);
    vector[N + 1] z = zeros_vector(N + 1);
    real sum_w = 0;
    for (ii in 1:N) {
      int i = N - ii + 1; 
      real n = i;
      real w = y[i] * inv_sqrt(n * (n + 1));
      sum_w += w;
      z[i] += sum_w;     
      z[i + 1] -= w * n;    
    }
    return z;
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
  int<lower=1> K_r; // must have at least one 
  int<lower=1> N_grp;
  array[K_r] int<lower=1,upper=N_grp> group;
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
  array[n_nonzero] int v = csr_extract_v(X_r_all);
  array[N_obs+N_cens+1] int u= csr_extract_u(X_r_all);

  // Random affect group sizes
  array[N_grp] int comp_sizes = component_size(group);

  if(K){
	// thin and scale the QR decomposition
	Q_ast = qr_thin_Q(append_row(X_obs, X_cens)) * sqrt((N_obs+N_cens) - 1);
	R_ast = qr_thin_R(append_row(X_obs, X_cens)) / sqrt((N_obs+N_cens) - 1);
	R_ast_inverse = inverse(R_ast);
  }
}
parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] thetas; // coefficients on Q_ast
  vector[(K_r - N_grp)] rand_betas_raw;
  vector<lower=0>[N_grp] rand_sigma;
  real<lower=0> sigma_ln_eDNA;
}

transformed parameters {
  vector[K] betas;
  vector[K_r] rand_betas;
  if(K){
	betas = R_ast_inverse * thetas; // coefficients on x
  }
  
  rand_betas = zero_sum_components(rand_betas_raw, comp_sizes);

  for(k in 1:K_r)
	rand_betas[k] = rand_betas[k] * rand_sigma[group[k]];

}

model {
  vector[N_obs + N_cens] mu_all =
	(K ? Q_ast * thetas : rep_vector(0.0, N_obs + N_cens)) +
	csr_matrix_times_vector(N_obs+N_cens, K_r, w, v, u, rand_betas) +
  /* (append_row(X_obs_r, X_cens_r) * rand_betas) + */
	(has_inter ? intercept[1] : 0.0);
  vector[N_obs] mu_obs = mu_all[1:N_obs];
  vector[N_cens ? N_cens : 1] mu_cens = N_cens ? mu_all[(N_obs+1):] : [0]';

  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);
  rand_sigma ~ gamma(2, 0.1);

  if(K) {
	for(k in 1:K)
	  betas[k] ~ normal(prior_mu[k], prior_sd[k]);
  }
  
  rand_betas_raw ~ std_normal();
  
  sigma_ln_eDNA ~ exponential(1);
  
  target += normal_lpdf(y_obs | mu_obs, sigma_ln_eDNA);
  if(N_cens)
	target += normal_lcdf(L | mu_cens, sigma_ln_eDNA);
}

generated quantities{
  vector[N_obs + N_cens] log_lik = rep_vector(0, N_obs + N_cens);
  {
    matrix[N_obs, K + K_r] X_obs_tmp = append_col(X_obs, X_obs_r);
    matrix[N_cens, K + K_r] X_cens_tmp = append_col(X_cens, X_cens_r);
    vector[K+K_r] betas_tmp = append_row(betas, rand_betas);
    if(N_obs) {
      for(n in 1:N_obs){
	log_lik[n] = normal_lpdf(y_obs[n] | (has_inter ? intercept[1] : 0.0) +
				 X_obs_tmp[n] * betas_tmp, sigma_ln_eDNA);
      }
    }
    
    if(N_cens) {
      for(n in 1:N_cens){
	log_lik[n+N_obs] = normal_lcdf(L[n] | (has_inter ? intercept[1] : 0) +
				     X_cens_tmp[n] * betas_tmp, sigma_ln_eDNA);
      }
    }
  }
}

