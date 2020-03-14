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
  
  int[] get_group_lengths(int[] groups){
	int n_grp = max(groups);
	int n = dims(groups)[1];
	int ans[n_grp] = rep_array(0, n_grp);
	for(i in 1:n){
	  ans[groups[i]] += 1;
	}
	return ans;
  }

  // Enforce sum to 0 on randoms
  vector raw_to_beta(vector raw_betas){
	int n = dims(raw_betas)[1];
	vector[n + 1] ans = append_row(raw_betas, -sum(raw_betas));

	return ans;		   
  }

  /*
	Not finished
	vector expand_raw_betas(vector raw_betas, int[] grps){
	int N = size(grps);
	int n_grps = max(grps);
	int group_lengths[n_grps] = get_group_lengths(grps);
	vector[N] betas;
	
	for(n in 1:N){
	  
	}
	
  }
  */
  
  // multiply the vector of random effects by the correct
  // variance param
  vector make_random_betas(vector raw_betas,
						   int[] grps,
						   vector random_sigmas){
	int n = dims(raw_betas)[1];
	vector[n] betas;
	
	for(i in 1:n)
	  betas[i] = raw_betas[n] * random_sigmas[grps[n]];
	
	return betas;  
  }
}

data{
  int N;
  int n_vars;
  matrix[N, n_vars] X;
  vector[N] y; // Cq 

  // stuff for random effects
  int<lower = 0, upper = 1> has_random;
  int n_rand; // number of columns of rand effects
  int n_grp;
  int groups[n_rand];
  matrix[N,n_rand] rand_x;

  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha;
  real std_curve_beta;
  real upper_Cq; // upper value that Cq can take
}

transformed data{
  // thin and scale the QR decomposition
  matrix[N, n_vars] Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  matrix[n_vars, n_vars] R_ast = qr_thin_R(X) / sqrt(N - 1);
  matrix[n_vars, n_vars] R_ast_inverse = inverse(R_ast);

  int group_lengths[has_random ? n_grp : 0];
  
  if(has_random)
	group_lengths = get_group_lengths(groups);
}

parameters{
  vector[n_vars] thetas;
  real<lower = 0> sigma_Cq;
  vector[has_random ? n_rand : 0] rand_betas_raw;
  vector<lower = 0>[has_random ? n_grp : 0] rand_sigma;
}
transformed parameters{
  vector[n_vars] betas = R_ast_inverse * thetas;
  vector[has_random ? n_rand : 0] rand_betas;
  if(has_random)
	rand_betas = make_random_betas(rand_betas_raw, groups, rand_sigma);
}

model{
  vector[N] ln_conc_hat = Q_ast * thetas;
  vector[N] Cq_hat;

  if(has_random)
	ln_conc_hat = ln_conc_hat + rand_x * rand_betas;
  
  // Priors
  thetas ~ normal(0 , 1);
  sigma_Cq ~ normal(0, 1);
  rand_sigma ~ normal(0, 1);
  
  // not sure if this works
  rand_betas_raw ~ normal(0, 1);
  
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
 
}
