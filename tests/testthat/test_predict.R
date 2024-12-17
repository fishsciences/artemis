context("Model predictions")

test_that("Model predictions",{
    dd = eDNA_lm(Cq ~ Distance_m, eDNA_data, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    ans = predict(dd)
    
    expect_is(ans, "eDNA_predict_lm")

    ans2 = predict(dd, interval = TRUE)
    expect_is(ans2, "eDNA_predict_lm")
    expect_true(ncol(ans2$ln_conc) == 2)
    expect_true(ncol(ans2$Cq_hat) == 2)
    expect_true(is.null(ans$Cq_star))

    ans3 = predict(dd, interval = TRUE, include_sigma = TRUE)
    expect_true(!is.null(ans3$Cq_star))

    x = data.frame(Distance_m = c(0,50))
 
    ans4 = predict(dd, newdata = x)

    
    
    expect_true(nrow(ans4$ln_conc) == 2)


    ans5 = predict(dd, include_sigma = TRUE)
    
})
