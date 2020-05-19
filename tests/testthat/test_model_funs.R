context("Model fitting")

ans = eDNA_lm(Cq ~ Distance, eDNA_data,
              std_curve_alpha = 21.2, std_curve_beta = -1.5)
ans2 = eDNA_lm(Cq ~ Distance + Volume, eDNA_data,
               std_curve_alpha = 21.2, std_curve_beta = -1.5)
ans_prior = eDNA_lm(Cq ~ Distance, eDNA_data,
              std_curve_alpha = 21.2, std_curve_beta = -1.5,
              prior_intercept = normal(-8,1),
              priors = normal(0, 1))

test_that("Fit the model with simple data",{
    expect_is(ans, "eDNA_model_lm")

    expect_true(all(slotNames(ans) %in% c("ln_conc", "Cq_star", "intercept", "betas",
                                          "sigma_Cq", "formula", "x", 
                                          "std_curve_alpha", "std_curve_beta",
                                          "upper_Cq", "stanfit")))

    aa = summary(ans)

    expect_is(aa, "data.frame")

    expect_is(ans, "eDNA_model_lm")

    expect_true(all(slotNames(ans) %in% c("ln_conc", "Cq_star", "intercept","betas",
                                          "sigma_Cq", "formula", "x", 
                                          "std_curve_alpha", "std_curve_beta",
                                          "upper_Cq", "stanfit")))

    aa = summary(ans)
})

test_that("Lmer", {
    ans2 = eDNA_lmer(Cq ~ Distance + Volume + (1|SampleID), eDNA_data,
                   std_curve_alpha = 21.2, std_curve_beta = -1.5, verbose = FALSE)

   
    summary(ans2)

    #Only run if multicore available
    if(FALSE)
        ans2 = eDNA_lmer(Cq ~ Distance + Volume + (1|SampleID),
                         eDNA_data,
                         std_curve_alpha = 21.2, std_curve_beta = -1.5, verbose = FALSE,
                         cores = floor(parallel::detectCores() / 2), iter = 2000)
    
})

test_that("lm with priors", {
    # This should still work
    # ans = eDNA_lm(Cq ~ Distance, eDNA_data,
    #               std_curve_alpha = 21.2, std_curve_beta = -1.5)

                  
    # d = eDNA_data
    # d$Distance = 5

    # Does not run - not sure why not
    # ans = eDNA_lm(Cq ~ Distance -1, d,
    #               std_curve_alpha = 21.2, std_curve_beta = -1.5,
    #               prior_intercept = normal(-15, 5),
    #               priors = normal(-2, 5))

})

test_that("Intercepts", {
    ans = eDNA_lm(Cq ~ 1,  eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(ans, "eDNA_model")

})

test_that("Varying measurement error", {
  ans = eDNA_lm(Cq ~ Distance + Volume, eDNA_data,
                std_curve_alpha = 21.2, std_curve_beta = -1.5,
                Cq_error_type = "varying")

  expect_is(ans, "eDNA_model")

  expect_error(eDNA_lm(Cq ~ Distance + Volume, eDNA_data,
                std_curve_alpha = 21.2, std_curve_beta = -1.5,
                Cq_error_type = "bob"))


})
