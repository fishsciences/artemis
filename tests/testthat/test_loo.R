context("LOO")

ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
              std_curve_alpha = 21.2, std_curve_beta = -1.5)

test_that("loo", {
    expect_is(loo(ans), "psis_loo")

})

test_that("waic", {
    expect_is(waic(ans), "waic")
})
