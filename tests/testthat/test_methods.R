context("Random effects")

test_that("Random effects", {
    ans2 = eDNA_lmer(Cq ~ Distance + Volume + (1|Distance), eDNA_data,
                     std_curve_alpha = 21.2, std_curve_beta = -1.5, verbose = FALSE)

    rr = ranef(ans2)

    expect_is(rr, "matrix")
    expect_true(nrow(rr) == 3)


})
