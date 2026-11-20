  vector[N_obs + N_cens] mu_all =
    (K ? Q_ast * thetas : rep_vector(0.0, N_obs + N_cens)) +
    (has_inter ? intercept[1] : 0.0);
