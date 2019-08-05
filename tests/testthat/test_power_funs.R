context("Power functions")

test_that("P-detect", {
    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 500),
                              betas = c(intercept = -10.6, Distance = -0.005, Volume = 0.01),
                              Cq_sd = 1, n_rep = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5))
    
    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 500, Volume = 50),
                 betas = c(intercept = -10.6, Distance = -0.005),
                 Cq_sd = 1, n_rep = 6, std_curve_alpha = 21.2, std_curve_beta = -1.5))

    res = est_p_detect(var_levels = c(intecept = 1, Distance = 300, Volume = 0),
                 betas = c(intercept = -10.6, Distance = -0.5, Volume = -0.1),
                 Cq_sd = 1, n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    ans = eDNA_lm(Cq ~ Distance, eDNA_samples,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)


    res = est_p_detect(var_levels = c(intecept = 1, Distance = 50),
                 betas = NULL, model_fit = ans, 
                 Cq_sd = 1, n_rep = 2, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == 1)

    res = est_p_detect(var_levels = c(intecept = 1, Distance = 50),
                 betas = NULL, model_fit = ans, 
                 Cq_sd = 1, n_rep = 2:6, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == length(2:6))


    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 300,
                                             Volume = c(0, 20)),
                              betas = c(intercept = -10.6, Distance = -0.5, Volume = -0.1),
                              Cq_sd = 1, n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5))
   
    
})

test_that("Power function", {
    
    ans = est_power(Cq ~ 1 + distance + (1|rep), vars,
                    betas = c(intercept = -10.6, distance = -0.05),
                    sigma_Cq = 1, sigma_rep = 0.5, std_curve_alpha = 21.2, std_curve_beta = -1.5)


})

