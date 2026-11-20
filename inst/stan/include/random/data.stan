  // Random effects
  int<lower=1> K_r; // must have at least one
  int<lower=1> N_grp;
  array[K_r] int<lower=1,upper=N_grp> group; // contiguous by group
  array[N_grp] int<lower=1> grp_size;
  real<lower=0> rand_sd; // rate of exponential prior on rand_sigma
