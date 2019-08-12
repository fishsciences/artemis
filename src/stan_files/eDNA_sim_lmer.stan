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
  
  int groups[n_rand];
  matrix[N,n_rand] rand_x;
  real upper_Cq; // upper value that Cq can take
  
  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha;
  real std_curve_beta;

  real<lower = 0> sigma_Cq; // sd on Cq - assumed to be only on measurement
  vector[n_vars] betas;
  vector<lower = 0>[n_grp] rand_sigma;

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
  vector[N] ln_conc;
  vector[N] Cq_star;
  vector[n_rand] rand_betas;

  for(i in 1:n_rand)
	rand_betas[i] = normal_rng(0, 1) * rand_sigma[groups[i]];

  ln_conc = X * betas + rand_x * rand_betas;
  
  for(n in 1:N){
	real Cq_hat = ln_conc[n] * std_curve_beta + std_curve_alpha;
	
	Cq_star[n] = normal_rng(Cq_hat, sigma_Cq);
	if(Cq_star[n] > upper_Cq)
	  Cq_star[n] = upper_Cq;
  }
}
