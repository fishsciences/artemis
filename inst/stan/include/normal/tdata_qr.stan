  // thin and scale the QR decomposition
  matrix[N_obs + N_cens, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  if (K) {
    Q_ast = qr_thin_Q(append_row(X_obs, X_cens)) * sqrt((N_obs + N_cens) - 1);
    R_ast = qr_thin_R(append_row(X_obs, X_cens)) / sqrt((N_obs + N_cens) - 1);
    R_ast_inverse = inverse(R_ast);
  }
