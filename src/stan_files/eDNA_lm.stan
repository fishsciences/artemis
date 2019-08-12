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

}

data{
  int N;
  int n_vars;

  matrix[N, n_vars] X;
  real upper_Cq; // upper value that Cq can take

  vector[N] y; // Cq 
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
  real<lower = 0> sigma_Cq;
}
transformed parameters{
  
}

model{
  vector[N] ln_conc_hat = Q_ast * thetas;
  vector[N] Cq_hat;
  
  // Priors
  thetas ~ normal(0 , 1);
  sigma_Cq ~ normal(0, 1);
  
  Cq_hat = ln_conc_hat * std_curve_beta + std_curve_alpha;
  
  for(n in 1:N){
	if(y[n] < upper_Cq) {
	  y[n] ~ normal(Cq_hat[n], sigma_Cq);  
	} else {
	  target += normal_lccdf(upper_Cq | Cq_hat[n], sigma_Cq);
	}
  }
}

generated quantities{
  vector[n_vars] betas = R_ast_inverse * thetas;

}
