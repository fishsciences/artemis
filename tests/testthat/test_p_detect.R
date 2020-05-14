context("P-detect")

test_that("Duplicate arg warnings", {
    expect_warning(dup_arg_warn("bob"), "bob")

    model_fit = eDNA_lm(Cq ~ Distance, eDNA_data,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_warning(est_p_detect(variable_levels = c(Distance = 1000),
                            model_fit = model_fit, 
                            Cq_sd = 1, 
                            std_curve_alpha = 21.2, std_curve_beta = -1.5,
                            n_rep = 12:30))

})

test_that("Multiple var levels", {
    expect_error(est_p_detect(variable_levels = data.frame(Intercept = 1, 
                                                           Distance = c(10, 1000)),
                              betas = c(Intercept = -12, Distance = -0.0001),
                              Cq_sd = 1, 
                              std_curve_alpha = 21.2, std_curve_beta = -1.5,
                              n_rep = 12:30))


    })
    
