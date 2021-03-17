context("P-detect")

test_that("Duplicate arg warnings", {
    expect_warning(dup_arg_warn("bob"), "bob")

    model_fit = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_warning(est_p_detect(variable_levels = c(Distance_m = 1000),
                            model_fit = model_fit, 
                            ln_eDNA_sd = 1, 
                            std_curve_alpha = 21.2, std_curve_beta = -1.5,
                            n_rep = 12:30))

})

test_that("Multiple var levels", {
    expect_error(est_p_detect(variable_levels = data.frame(Intercept = 1, 
                                                           Distance_m = c(10, 1000)),
                              betas = c(Intercept = -12, Distance_m = -0.0001),
                              ln_eDNA_sd = 1, 
                              std_curve_alpha = 21.2, std_curve_beta = -1.5,
                              n_rep = 12:30))


    })
    
test_that("Zero-inflated", {
    model_fit = eDNA_lm(Cq ~ scale(Distance_m), eDNA_data,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5)

    est_p_detect(variable_levels = c(Distance_m = 2),
                 model_fit = model_fit, 
                 std_curve_alpha = 21.2, std_curve_beta = -1.5,
                 n_rep = 1:3)
    est_p_detect(variable_levels = c(Distance_m = 0),
                 model_fit = model_fit, 
                 std_curve_alpha = 21.2, std_curve_beta = -1.5,
                 n_rep = 1:3)

})
