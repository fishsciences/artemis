context("User Priors")

test_that("User priors: lm", {
    # Very basic test - need more thorough tests
    # no priors
    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(ans, "eDNA_model_lm")
    # Priors
    ans2 = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                   std_curve_alpha = 21.2, std_curve_beta = -1.5,
                   prior_intercept = normal(-12, 10),
                   priors = normal(0,1))
    
    expect_is(ans2, "eDNA_model_lm")
    
})

test_that("User priors: random effect SD", {
    ml = gen_model_list_lmer(Cq ~ Distance_m + (1|FilterID), eDNA_data)
    md = prep_data.lmer(ml, 21.2, -1.5, Cq_upper = 40, prior_int = normal(),
                        prior_b = normal(), rand_sd = exponential(10))
    expect_equal(md$rand_sd, 10)
    expect_error(prep_data.lmer(ml, 21.2, -1.5, Cq_upper = 40, prior_int = normal(),
                                prior_b = normal(), rand_sd = normal()))

    # a tight prior should shrink rand_sigma relative to a wide one
    tight = eDNA_lmer(Cq ~ Distance_m + (1|FilterID), eDNA_data,
                      std_curve_alpha = 21.2, std_curve_beta = -1.5,
                      prior_random_variance = exponential(10), seed = 1234)
    wide = eDNA_lmer(Cq ~ Distance_m + (1|FilterID), eDNA_data,
                     std_curve_alpha = 21.2, std_curve_beta = -1.5,
                     prior_random_variance = exponential(0.1), seed = 1234)
    expect_true(median(tight@fit$draws("rand_sigma")) <
                median(wide@fit$draws("rand_sigma")))
})
