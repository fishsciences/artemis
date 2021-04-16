/*
  Stan model for simulating data from an eDNA process

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

}

data{
  int N;
  int n_vars;
  int n_rand; // number of columns of rand effects
  int n_grp;
  matrix[N, n_vars] X;
  
  int<lower = 0, upper = 1> has_random;
  int groups[n_rand];
  matrix[N,n_rand] rand_x;
  vector<lower = 0>[n_grp] rand_sigma;
  
  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha[N];
  real std_curve_beta[N];
  real upper_Cq; // upper value that Cq can take

  // zero-inflated
  real<lower = 0, upper = 1> p_zero;
  
  real<lower = 0> sigma_ln_eDNA; // sd on Cq - assumed to be only on measurement
  vector[n_vars] betas;

}

transformed data{
}

parameters{
}

transformed parameters{ 
}

model{
}

generated quantities{
  vector[N] Cq_star;
  vector[N] ln_conc_star;
  vector[N] ln_conc_hat;
  vector[has_random ? n_rand : 0] rand_betas;

  ln_conc_star = X * betas;

  if(has_random){
	for(i in 1:n_rand)
	  rand_betas[i] = normal_rng(0, 1) * rand_sigma[groups[i]];
	
	ln_conc_star = ln_conc_star + rand_x * rand_betas;
  }
  
  for(n in 1:N){
	ln_conc_hat[n] = normal_rng(ln_conc_star[n], sigma_ln_eDNA);

	Cq_star[n] = ln_conc_hat[n] * std_curve_beta[n] + std_curve_alpha[n];
	if(Cq_star[n] > upper_Cq)
	  Cq_star[n] = upper_Cq;

	// zero-inflated part
	if(bernoulli_rng(p_zero))
	  Cq_star[n] = upper_Cq;

  }
}
