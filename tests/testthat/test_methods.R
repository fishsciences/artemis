context("Random effects")

test_that("Random effects", {
    ans2 = eDNA_lmer(Cq ~ Distance_m + Volume_mL + (1|Distance_m), eDNA_data,
                     std_curve_alpha = 21.2, std_curve_beta = -1.5)

    rr = ranef(ans2)
    n_lev = length(unique(eDNA_data$Distance_m))

    expect_is(rr, "matrix")
    expect_true(nrow(rr) == n_lev)
    expect_true(ncol(rr) == 3)
    expect_true(rownames(rr)[1] == "(Intercept):Distance_m[10]")

    rr = ranef(ans2, FUN = mean)
    expect_true(all(dim(rr) == c(n_lev, 1)))

    res = summary(ans2)
    expect_true(tail(res$variable, 1) == "sd((Intercept):Distance_m)")
})

test_that("Random effects: zero-inflated and count models", {
    d = eDNA_data
    d$count = ceiling(12.9 - 0.367 * d$Cq)
    d$count[d$Cq == 40] = 0
    d$count[d$count < 0] = 1

    m_zip = eDNA_zinf_lmer(Cq ~ Distance_m + (1|Distance_m) | 1, d,
                           std_curve_alpha = 21.2, std_curve_beta = -1.5)
    m_count = eDNA_count_lmer(count ~ 1 + (1|Distance_m), d)

    for(m in list(m_zip, m_count)){
        rr = ranef(m)
        expect_true(nrow(rr) == length(unique(d$Distance_m)))
        expect_true(all(grepl("^\\(Intercept\\):Distance_m\\[", rownames(rr))))
        expect_true("sd((Intercept):Distance_m)" %in% summary(m)$variable)
    }
})

