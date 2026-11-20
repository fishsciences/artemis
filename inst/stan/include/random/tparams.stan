  vector[K_r] rand_betas =
    zero_sum_components(rand_betas_raw, grp_size) .* rand_sigma[group];
