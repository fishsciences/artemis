/* Stan model for a Poisson distributed linear model
   M. Espe
*/


data {
  int<lower=0> N;
  int<lower=0> K;

  matrix[N,K] X;

  // prior parameters - user provided
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;  

  // for intercept-less models
  int<lower=0,upper=1> has_inter;

  array[N] int<lower=0> y;
}

parameters {
  array[has_inter ? 1 : 0] real intercept;
  vector[K] betas;      // coefficients on Q_ast
}

transformed parameters {

}

model {
  // priors
  intercept ~ normal(prior_int_mu, prior_int_sd);

  if(K) {
    for(k in 1:K)
      betas[k] ~ normal(prior_mu[k], prior_sd[k]);
  }
  
  // Likelihood
  y ~ poisson_log_glm(X, (has_inter ? intercept[1] : 0.0), betas);
}

generated quantities {
  vector[N] log_lik;

  for(n in 1:N){
	log_lik[n] = poisson_log_lpmf(y[n] | X[n] * betas + (has_inter ? intercept[1] : 0.0));
  }
    
}
