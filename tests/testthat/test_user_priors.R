context("User Priors")

test_that("User priors: lm", {
    # Very basic test - need more thorough tests
    # no priors
    ans = eDNA_lm(Cq ~ Distance, eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(ans, "eDNA_model_lm")
    # Priors
    ans2 = eDNA_lm(Cq ~ Distance, eDNA_data,
                   std_curve_alpha = 21.2, std_curve_beta = -1.5,
                   betas_prior_mu = c(0,0), betas_prior_sd = c(1,1))
    expect_is(ans2, "eDNA_model_lm")
    
})
