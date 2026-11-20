  // zero-inflation components
  int<lower=0> K_z; // predictors
  matrix[N_cens, K_z] X_z;
  matrix[N_obs, K_z] X_nz;
  int<lower=0,upper=1> has_nz_inter;
