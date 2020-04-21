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

  int[] get_group_end(int[] groups){
	int n_grp = max(groups);
	int n = dims(groups)[1];
	int ans[n_grp] = rep_array(0, n_grp);
	int i = 1;
	int pos = 1;
	
	while(pos < n_grp){
	  if(groups[i + 1] != groups[i]){
		ans[pos] = i;
		pos += 1;
	  }
	  i += 1;
	}
	ans[n_grp] = n; // the last group ends at the end of the vector

	return ans;
  }

  vector groups_sum_to_zero(vector raw_betas, int[] group_ends){
	
	int n_grps = size(group_ends);
	int start = 1;
	vector[dims(raw_betas)[1]] ans = raw_betas;

	for(i in 1:n_grps){
	  ans[group_ends[i]] = -sum(ans[start:(group_ends[i] - 1)]);
	  start = group_ends[i] + 1;
	}
	
	return ans;	
  }

  // Enforce sum to 0 on randoms
  vector raw_to_beta(vector raw_betas){
	int n = dims(raw_betas)[1];
	vector[n + 1] ans = append_row(raw_betas, -sum(raw_betas));

	return ans;		   
  }
  
  // multiply the vector of random effects by the correct
  // variance param
  vector make_random_betas(vector raw_betas,
						   int[] grps,
						   vector random_sigmas){
	int n = dims(raw_betas)[1];
	vector[n] betas;
	
	for(i in 1:n)
	  betas[i] = raw_betas[i] * random_sigmas[grps[i]];
		
	return betas;  
  }

  
}

data{
  int N;
  int n_vars;
  matrix[N, n_vars] X;
  vector[N] y; // Cq 

  // intercept stuff
  int<lower = 0, upper = 1> has_intercept;
  real intercept_mu;
  real intercept_sd;
  int<lower = 0, upper = 1> hs_prior; // 0 = normal, 1 = hs_plus

  // stuff for random effects
  int<lower = 0, upper = 1> has_random;
  int n_rand; // number of columns of rand effects
  int n_grp;
  int groups[n_rand];
  matrix[N,n_rand] rand_x;

  // user priors
  int<lower = 0, upper = 1> has_prior;
  vector[n_vars] prior_mu;
  vector[n_vars] prior_sd;  
  
  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha;
  real std_curve_beta;
  real upper_Cq; // upper value that Cq can take

  // for HS prior - from Piironen and Vehtari 2018
  real<lower = 0> global_scale; // scale for the half -t prior for tau
  real<lower = 1> nu_global; // degrees of freedom for the half -t priors for tau
  real<lower = 1> nu_local; // degrees of freedom for the half - t priors for lambdas
  real<lower = 0> slab_scale; // slab scale for the regularized horseshoe
  real<lower = 0> slab_df; // slab degrees of freedom for the regularized horseshoe
}

transformed data{
  // thin and scale the QR decomposition
  matrix[N, n_vars] Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  matrix[n_vars, n_vars] R_ast = qr_thin_R(X) / sqrt(N - 1);
  matrix[n_vars, n_vars] R_ast_inverse = inverse(R_ast);

  int group_ends[has_random ? n_grp : 0];
  
  if(has_random){
	group_ends = get_group_end(groups);
  }
}

parameters{
  real alpha[has_intercept ? 1 : 0];
  vector[n_vars] thetas;
  real<lower = 0> sigma_Cq;
  vector[has_random ? n_rand : 0] rand_betas_raw;
  vector<lower = 0>[has_random ? n_grp : 0] rand_sigma;
  real<lower=0> tau[hs_prior ? 1 : 0]; // global shrinkage parameter
  vector<lower=0>[hs_prior ? n_vars : 0] lambda; // local shrinkage parameter
  real<lower=0> caux[hs_prior ? 1 : 0];
}

transformed parameters{
  vector[n_vars] betas;
  vector[has_random ? n_rand : 0] rand_betas;
  real<lower=0> c[hs_prior ? 1 : 0];

  vector<lower=0>[hs_prior ? n_vars : 0] lambda_tilde;

  if(hs_prior){
	c[1] = slab_scale * sqrt(caux[1]); // slab scale
	lambda_tilde = sqrt(c[1] ^ 2 * square(lambda) ./
						(c[1] ^ 2 + tau[1] ^ 2 * square(lambda)));

	betas = R_ast_inverse * thetas .* lambda_tilde * tau[1];
  } else {
	betas = R_ast_inverse * thetas;
  }
  
  if(has_random){
	rand_betas = groups_sum_to_zero(rand_betas_raw, group_ends);
  	rand_betas = make_random_betas(rand_betas, groups, rand_sigma);
  }
}

model{
  vector[N] ln_conc_hat = has_intercept ? alpha[1] + X * betas :
	                                      X * betas;
  vector[N] Cq_hat;

  if(has_random)
	ln_conc_hat = ln_conc_hat + rand_x * rand_betas;
  
  // Priors
  if(has_intercept)
	alpha ~ normal(intercept_mu, intercept_sd);
  
  sigma_Cq ~ normal(0, 1);
  rand_sigma ~ normal(0, 1);

  // HS+
  if(hs_prior) {
	lambda ~ student_t(nu_local, 0, 1);
	tau ~ student_t(nu_global, 0, global_scale * 1); // was sigma
	caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  }
  
  // not sure if this works
  rand_betas_raw ~ normal(0, 1);

  if(has_prior){
	betas ~ normal(prior_mu, prior_sd);
  } else {
	thetas ~ normal(0, 1);
  }
  
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
