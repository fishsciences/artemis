context("Model fitting")

test_that("Fit the model with simple data",{
    ans = eDNA_lm(Cq ~ Distance, eDNA_samples, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_is(ans, "eDNA_model_lm")

    expect_true(all(slotNames(ans) %in% c("ln_conc", "Cq_star", "betas",
                                          "sigma_Cq", "formula", "x", 
                                          "std_curve_alpha", "std_curve_beta",
                                          "upper_Cq", "stanfit")))

    aa = summary(ans)

    expect_is(aa, "data.frame")
    
    
})
