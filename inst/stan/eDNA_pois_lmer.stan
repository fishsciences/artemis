/* Stan model for a Poisson distributed linear model
   M. Espe
*/

functions {
  // Adapted from https://github.com/mitzimorris/geomed_2024/blob/main/stan/bym2_multicomp.stan
  /**
   * Component-wise constrain sum-to-zero vectors
   *
   * @param phi unconstrained vector of zero-sum slices
   * @param sizes component sizes; components must be contiguous
   * @return vector phi_ozs, the vector whose slices sum to zero
   */
  vector zero_sum_components(vector phi, array[] int sizes) {
    vector[sum(sizes)] phi_ozs;
    int idx_phi = 1;
    int idx_ozs = 1;
    for (i in 1:size(sizes)) {
      int n = sizes[i];
      phi_ozs[idx_ozs : (idx_ozs + n - 1)] = 
        sum_to_zero_constrain(segment(phi, idx_phi, n - 1));
      idx_phi += n - 1;
      idx_ozs += n;
    }
    return phi_ozs;
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;

  matrix[N,K] X;

  // prior parameters - user provided
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;
  real<lower=0> rand_sd; // rate of exponential prior on rand_sigma

  // for intercept-less models
  int<lower=0,upper=1> has_inter;

  // Random effects
  int<lower=1> K_r; // must have at least one 
  int<lower=1> N_grp;
  array[K_r] int<lower=1,upper=N_grp> group; // contiguous by group
  array[N_grp] int<lower=1> grp_size;
  matrix[N, K_r] X_r;
    
  array[N] int<lower=0> y;
  
}

parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] betas;      // coefficients on Q_ast
  vector[(K_r - N_grp)] rand_betas_raw;
  vector<lower=0>[N_grp] rand_sigma;
}

transformed parameters {
  vector[K_r] rand_betas =
	zero_sum_components(rand_betas_raw, grp_size) .* rand_sigma[group];
}

model {
  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);
  rand_sigma ~ exponential(rand_sd);
  rand_betas_raw ~ std_normal();
  
  if(K) {
    for(k in 1:K)
      betas[k] ~ normal(prior_mu[k], prior_sd[k]);
  }
  
  // Likelihood
  y ~ poisson_log_glm(append_col(X, X_r),
					  (has_inter ? intercept[1] : 0.0),
					  append_row(betas, rand_betas));
}

generated quantities {
  vector[N] log_lik;

  for(n in 1:N){
	log_lik[n] = poisson_log_lpmf(y[n] |
								  append_col(X[n], X_r[n]) *
								  append_row(betas, rand_betas) +
								  (has_inter ? intercept[1] : 0.0));
  }
    
}
