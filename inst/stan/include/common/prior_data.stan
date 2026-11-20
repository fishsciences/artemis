  // prior parameters - user provided (requires K)
  real prior_int_mu;
  real<lower=0> prior_int_sd;
  vector[K] prior_mu;
  vector<lower=0>[K] prior_sd;

  // for intercept-less models
  int<lower=0,upper=1> has_inter;
