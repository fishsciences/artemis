/*
  Stan model for conducting a power analysis 
  for the eDNA project

  M. Espe
  March 2019

  The model is a measurement error model, where the 
  observed response is Cq, which a measurement of eDNA 
  concentrations + error.

  Since Stan does not allow the generation of random data 
  outside of the generated quantities block, we need to generate data
  in R, and then feed it to Stan to get estimates and calculate the 
  power for the test.
 */

functions{
  // Calculate standard curve conversion
  vector ln_std_curve(vector conc, real std_curve_alpha, real std_curve_beta){
	return std_curve_beta * log(conc) + std_curve_alpha;
  }
}

data{
  int N;
  int n_vars;
  int n_rand_var; // number of columns of rand effects
  int n_rand_total; // total number of random effects
  int<upper = n_rand_var> rand_var_shared[n_rand_total]; // idx to map which are shared
  matrix[N, n_vars] X;
  int rand_id[N, n_rand_var];
  vector[N] Cq;
  real upper_Cq; // upper value that Cq can take
  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha;
  real std_curve_beta;
}

transformed data{
  // QR Decomp - from Stan User Guide
  matrix[n_vars, n_vars] R = qr_thin_R(X);
  real s = sqrt(N - 1.0);
  matrix[n_vars, n_vars] R_ast = R / s;
  matrix[N, n_vars] Q_ast = qr_thin_Q(X) * s;
  matrix[n_vars, n_vars] R_ast_inverse = inverse(R_ast);
}

parameters{
  vector[n_vars] thetas;
  vector[n_rand_total] rand_betas_raw;
  real<lower = 0> sigma_Cq;
  vector<lower = 0>[n_rand_var] rand_sigma;
}
transformed parameters{
  
}

model{
  vector[N] ln_conc_hat = Q_ast * thetas;
  vector[N] Cq_hat;
  
  
  for(i in 1:n_rand_var)
	for(n in 1:N)
	  ln_conc_hat[n] += rand_betas_raw[rand_id[n,i]] * rand_sigma[i];

  // Priors
  thetas ~ normal(0 , 1);
  rand_betas_raw ~ normal(0, 1);
  sigma_Cq ~ normal(0, 1);
  rand_sigma ~ normal(0, .1);
  
  Cq_hat = ln_std_curve(exp(ln_conc_hat), std_curve_alpha, std_curve_beta);
  
  for(n in 1:N){
	if(Cq[n] < upper_Cq) {
	  Cq[n] ~ normal(Cq_hat[n], sigma_Cq);  
	} else {
	  target += normal_lccdf(upper_Cq | Cq_hat[n], sigma_Cq);
	}
  }
}

generated quantities{
  vector[n_vars] betas = R_ast_inverse * thetas;
  /* vector */
}
