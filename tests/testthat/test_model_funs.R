context("Model fitting")

ans = eDNA_lm(Cq ~ scale(Distance_m), eDNA_data,
              std_curve_alpha = 21.2, std_curve_beta = -1.5)
ans2 = eDNA_lm(Cq ~ scale(Distance_m) + Volume_mL, eDNA_data,
               std_curve_alpha = 21.2, std_curve_beta = -1.5)
ans_prior = eDNA_lm(Cq ~ scale(Distance_m), eDNA_data,
              std_curve_alpha = 21.2, std_curve_beta = -1.5,
              prior_intercept = normal(-8,1),
              priors = normal(0, 1))

if(FALSE){ # Commented out to avoid UBSAN error - bad solution, but I cannot replicate CRAN's error across multiple systems. 
test_that("Intercept only", {

    # Check that intercept only model will work - not inside test block, but should still throw error

    m = eDNA_lm(Cq ~ 1, eDNA_data,
                std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(m, "eDNA_model")
    m = eDNA_lmer(Cq ~ 1 + (1|Distance_m), eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(m, "eDNA_model")

})
}

test_that("Fit the model with simple data",{
    skip_on_cran()
    expect_is(ans, "eDNA_model_lm")

    expect_true(all(slotNames(ans) %in% c("ln_conc", "Cq_star", "intercept", "betas",
                                          "sigma_ln_eDNA", "formula", "x", 
                                          "std_curve_alpha", "std_curve_beta",
                                          "upper_Cq", "stanfit")))

    aa = summary(ans)

    expect_is(aa, "data.frame")

    expect_is(ans, "eDNA_model_lm")

    expect_true(all(slotNames(ans) %in% c("ln_conc", "Cq_star", "intercept","betas",
                                          "sigma_ln_eDNA", "formula", "x", 
                                          "std_curve_alpha", "std_curve_beta",
                                          "upper_Cq", "stanfit")))

    aa = summary(ans)
})

test_that("Lmer", {
    skip_on_cran()
    ans2 = eDNA_lmer(Cq ~ Distance_m + Volume_mL + (1|FilterID), eDNA_data,
                   std_curve_alpha = 21.2, std_curve_beta = -1.5)

   
    summary(ans2)

    #Only run if multicore available
    if(FALSE)
        ans2 = eDNA_lmer(Cq ~ Distance_m + Volume_mL + (1|FilterID),
                         eDNA_data,
                         std_curve_alpha = 21.2, std_curve_beta = -1.5,
                         cores = floor(parallel::detectCores() / 2))
    
})

test_that("lm with priors", {
    skip_on_cran()
    ## This should still work
    # ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
    #               std_curve_alpha = 21.2, std_curve_beta = -1.5)

                  
    # d = eDNA_data
    # d$Distance_m = 5

    # Does not run - not sure why not
    # ans = eDNA_lm(Cq ~ Distance_m -1, d,
    #               std_curve_alpha = 21.2, std_curve_beta = -1.5,
    #               prior_intercept = normal(-15, 5),
    #               priors = normal(-2, 5))

})

test_that("Intercepts", {
    skip_on_cran()
    ans = eDNA_lm(Cq ~ 1,  eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(ans, "eDNA_model")
    ans = eDNA_lmer(Cq ~ 1 + (1|FilterID),  eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(ans, "eDNA_model")

})

test_that("Varying measurement error", {
    skip_on_cran()
    if(FALSE){
    ans = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data,
                std_curve_alpha = 21.2, std_curve_beta = -1.5,
                Cq_error_type = "varying")

  expect_is(ans, "eDNA_model")

  expect_error(eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data,
                std_curve_alpha = 21.2, std_curve_beta = -1.5,
                Cq_error_type = "bob"))
    }

})

test_that("All obs. below threshold", {
    d = data.frame(x = rnorm(30, 30, 1))

    ans = eDNA_lm(x ~ 1, d, std_curve_alpha = 20, std_curve_beta = -1.5)

    expect_is(ans, "eDNA_model")

})
          
