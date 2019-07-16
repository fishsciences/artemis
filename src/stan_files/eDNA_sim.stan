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
  int<lower=0, upper=1> has_rand;
  int n_rand_var; // number of columns of rand effects
  int n_rand_total; // total number of random effects
  int<upper = n_rand_var> rand_var_shared[n_rand_total]; // idx to map which are shared
  int rand_id[N, n_rand_var];

  matrix[N, n_vars] X;
  int groups[ has_rand ? N : 0, n_rand_var];
  real upper_Cq; // upper value that Cq can take

  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha;
  real std_curve_beta;
  real<lower = 0> sigma_Cq; // sd on Cq - assumed to be only on measurement
  vector[n_vars] betas;
  vector<lower = 0>[n_rand_var] rand_sigma;
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
  /* vector */
  vector[N] ln_conc_fixed = X * betas;
  vector[N] Cq_star;

  if(has_rand) {
	real rand_beta_raw[n_rand_total] = normal_rng(rep_vector(0, n_rand_total), 1);

	for(i in 1:n_rand_var)
	  for(n in 1:N)
		ln_conc_fixed[n] += rand_beta_raw[rand_id[n,i]] * rand_sigma[i];
  }
  
  for(n in 1:N){
	real Cq_hat = ln_conc_fixed[n] * std_curve_beta + std_curve_alpha;

	Cq_star[n] = normal_rng(Cq_hat, sigma_Cq);
	if(Cq_star[n] > upper_Cq)
	  Cq_star[n] = upper_Cq;
  }
}
