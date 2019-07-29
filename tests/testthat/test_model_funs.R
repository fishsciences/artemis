context("Model fitting")

test_that("Fit the model with simple data",{
    ans = eDNA_lm(Cq ~ Distance, eDNA_samples, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(ans, "eDNA_model_lm")
 
    ans = eDNA_lm(Cq ~ Distance, eDNA_samples, fit_fun = rstan::optimizing,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
})
