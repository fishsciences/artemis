context("Multiple Std. Curves")

test_that("Multi-curves: model", {
    # TODO: add numerical checks
    n = nrow(eDNA_data)
    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = rep(21.2,n) , std_curve_beta = rep(-1.5, n))

    expect_is(ans, "eDNA_model")

    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = rnorm(n, 21.2, 0.025) ,
                  std_curve_beta = rnorm(n, -1.5, 0.025))
    
    expect_is(ans, "eDNA_model")
})

test_that("Multi-curves: methods", {
    n = nrow(eDNA_data)

    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = rnorm(n, 22.2, 0.025) ,
                  std_curve_beta = rnorm(n, -1.5, 0.025))

    res = predict(ans)
    # with random curves, we should not see any models
    expect_true(!any(duplicated(res)))
    
    expect_error(est_p_detect(10, model_fit = ans))

    ans2 = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = rep(21.2,n) , std_curve_beta = rep(-1.5, n))

    res2 = est_p_detect(10, model_fit = ans2)

    expect_is(res2, "eDNA_p_detect")   

})


