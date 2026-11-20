  // random effect X to sparse matrix for efficiency
  matrix[N_obs + N_cens, K_r] X_r_all = append_row(X_obs_r, X_cens_r);
  int n_nonzero = num_nonzero(X_r_all);
  vector[n_nonzero] w = csr_extract_w(X_r_all);
  array[n_nonzero] int v = csr_extract_v(X_r_all);
  array[N_obs + N_cens + 1] int u = csr_extract_u(X_r_all);
