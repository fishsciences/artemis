context("Model predictions")

test_that("Model predictions",{
    dd = eDNA_lm(Cq ~ Distance, eDNA_samples, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    ans = predict(dd)
    
    expect_is(ans, "eDNA_predict_lm")
 
    ans = eDNA_lm(Cq ~ Distance, eDNA_samples, fit_fun = rstan::optimizing,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
})
