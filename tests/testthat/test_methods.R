context("Random effects")

test_that("Random effects", {
    ans2 = eDNA_lmer(Cq ~ Distance_m + Volume_mL + (1|Distance_m), eDNA_data,
                     std_curve_alpha = 21.2, std_curve_beta = -1.5)

    rr = ranef(ans2)

    expect_is(rr, "matrix")
    expect_true(nrow(rr) == length(unique(eDNA_data$Distance_m)))


})
