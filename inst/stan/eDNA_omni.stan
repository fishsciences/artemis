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

  The data generating process is assumed to be:

  Cq ~ mix( normal(Cq_hat, sigma) [, Upper_Cq] , 0 ) // truncated
  0 ~ bernoulli(p_zero)
  Cq_hat = std_curve_beta * ln_conc_eDNA + std_curve_alpha
  ln_conc_eDNA = alpha + beta * X + U
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
  int n_vars; // excludes intercept
  matrix[N, n_vars] X;
  vector[N] y; // Cq 

  // stuff for random effects
  int<lower = 0, upper = 1> has_random;
  int n_rand; // number of columns of rand effects
  int n_grp;
  int groups[n_rand];
  matrix[N,n_rand] rand_x;

  // intercepts
  int<lower = 0, upper = 1> has_inter;
  real prior_int_mu;
  real<lower = 0> prior_int_sd;
  
  // user priors vs. QR
  int<lower = 0, upper = 1> use_qr;
  vector[n_vars] prior_mu;
  vector[n_vars] prior_sd;  

  // alpha and beta of the ln_conc -> Cq conversion according to
  // beta * log(conc) + alpha
  real std_curve_alpha[N];
  real std_curve_beta[N];
  real upper_Cq; // upper value that Cq can take

  // measurement error model
  int<lower = 0, upper = 1> sd_vary; // 0 = fixed, 1 = varying
  real<lower = 0, upper = upper_Cq> center_Cq; // CQ for centering 

  // For zero-inflation
  real<lower = 0, upper = 1> p_zero;

  // for vectorized sampling
  int n_below;
}

transformed data{
  matrix[N, n_vars] Xc; // centered X
  vector[n_vars] mean_x; // x means

  // thin and scale the QR decomposition
  matrix[N, n_vars] Q_ast; 
  matrix[n_vars, n_vars] R_ast; 
  matrix[n_vars, n_vars] R_ast_inverse; 

  int group_ends[has_random ? n_grp : 0];

  int n_above = N - n_below;
    
  if(has_random){
	group_ends = get_group_end(groups);
  }

  if(n_vars > 0) {
	for(i in 1:n_vars){ // either start at 1 (no inter) or 2 (inter)
	  mean_x[i] = mean(X[,i]);
	  Xc[,i] = X[,i] - mean_x[i];
	}
	
    Q_ast = qr_thin_Q(Xc) * sqrt(N - 1);      
	R_ast = qr_thin_R(Xc) / sqrt(N - 1);
	R_ast_inverse = inverse(R_ast);
  }
}

parameters{
  vector[has_inter ? 1 : 0] temp_intercept;
  vector[n_vars] thetas;
  real sigma_Cq_raw;
  vector[has_random ? n_rand : 0] rand_betas_raw;
  vector<lower = 0>[has_random ? n_grp : 0] rand_sigma;
  vector[sd_vary ? 1 : 0] sd_slope_location;
  vector<lower = 0>[sd_vary ? 1 : 0] sd_slope_scale;
  /* real<lower = 0, upper = 1> p_zero; */
}

transformed parameters{
  vector[n_vars] betas;
  vector[has_random ? n_rand : 0] rand_betas;
  
  if(has_random){
	rand_betas = groups_sum_to_zero(rand_betas_raw, group_ends);
  	rand_betas = make_random_betas(rand_betas, groups, rand_sigma);
  }

  if(n_vars > 0)
	betas = R_ast_inverse * thetas;
}

model{
  vector[N] ln_conc_hat;
  vector[N] Cq_hat;
  vector[N] sigmas = rep_vector(sigma_Cq_raw, N);
  real sd_slope_temp = 0;
  
  if(n_vars > 0) {
	ln_conc_hat = Q_ast * thetas + (has_inter ? temp_intercept[1] : 0.0);
  } else {
	ln_conc_hat = rep_vector(temp_intercept[1], N);
  }

  if(has_random)
	ln_conc_hat = ln_conc_hat + rand_x * rand_betas;
  
  for(n in 1:N)
	Cq_hat[n] = ln_conc_hat[n] * std_curve_beta[n] + std_curve_alpha[n];
  
  // Priors
  sigma_Cq_raw ~ normal(0.0, 0.5);
  rand_sigma ~ exponential(2);

  /* p_zero ~ normal(0.08, 0.01); */
  
  // not sure if this works
  rand_betas_raw ~ std_normal();

  if(has_inter)
	temp_intercept ~ normal(prior_int_mu, prior_int_sd);

  // For variable measurement error by Cq_hat
  if(sd_vary){
	sigmas = exp(sigmas + ((Cq_hat - center_Cq) * sd_slope_location[1] * sd_slope_scale[1]));
	sd_slope_location ~ std_normal();
	// 30 to 40 CQ = +1 sd - seems like a reasonable prior
	sd_slope_scale ~ normal(0.0, 0.1);
  }
    
  if(n_vars > 0) {
	if(use_qr){ 
	  for(i in 1:n_vars)
		betas[i] ~ normal(prior_mu[i], prior_sd[i]);
	} else {
	  thetas ~ std_normal();
	}
  }

  // convert to natural scale
  sigmas = exp(sigmas);
  
  if(n_below)
	target += normal_lpdf(head(y, n_below) |
						  head(Cq_hat, n_below), head(sigmas, n_below)) +
	  bernoulli_lpmf(0 | p_zero) * n_below;   
  if(n_above)
	target += log_sum_exp(bernoulli_lpmf(1 | p_zero) * n_above,
						  bernoulli_lpmf(0 | p_zero) * n_above +
						  normal_lccdf(rep_vector(upper_Cq, n_above) |
									   tail(Cq_hat, n_above),
									   tail(sigmas, n_above)));
  
}

generated quantities{
  real intercept[has_inter ? 1 : 0];
  real sd_slope[sd_vary ? 1 : 0];
  vector[N] log_lik;
  real sigma_Cq = exp(sigma_Cq_raw);
  
  if(n_vars > 0){
	intercept[1] = temp_intercept[1] - dot_product(mean_x, betas);
  } else {
 	intercept[1] = temp_intercept[1];
  }
  if(sd_vary)
	sd_slope[1] = sd_slope_location[1] * sd_slope_scale[1];

  // copied from above - needed to avoid saving all this as parameters
  // slightly inefficient, but saves on processing time later
  {
	  vector[N] ln_conc_hat_tmp;
	  vector[N] Cq_hat_tmp;
	  real sd_slope_temp_tmp = 0;

	  if(n_vars > 0) {
		ln_conc_hat_tmp = Q_ast * thetas + (has_inter ? temp_intercept[1] : 0.0);
	  } else {
		ln_conc_hat_tmp = rep_vector(temp_intercept[1], N);
	  }
	  
	  if(has_random)
		ln_conc_hat_tmp = ln_conc_hat_tmp + rand_x * rand_betas;
  
	  sd_slope_temp_tmp = sd_vary ? sd_slope_location[1] * sd_slope_scale[1] : 0.0;

	  for(n in 1:N)
		Cq_hat_tmp[n] = ln_conc_hat_tmp[n] * std_curve_beta[n] + std_curve_alpha[n];
	  
	  // For WAIC
	  for(n in 1:N){
		if(y[n] < upper_Cq) {
		  log_lik[n] = normal_lpdf(y[n] |
								   Cq_hat_tmp[n],
								   exp(sigma_Cq_raw +
									   (Cq_hat_tmp[n] - center_Cq) * sd_slope_temp_tmp)) +
			bernoulli_lpmf(0 | p_zero);   
		} else {
		  log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p_zero),
								   bernoulli_lpmf(0 | p_zero) +
								   normal_lccdf(upper_Cq | Cq_hat_tmp[n],
												exp(sigma_Cq_raw +
													((upper_Cq - center_Cq) * sd_slope_temp_tmp))));
		}
	  }
  }
}
