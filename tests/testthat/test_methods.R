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

test_that("Random slopes: effects sum to zero within each term", {
    skip_on_cran()
    d = eDNA_data
    # centered, so the slope column has negative values
    d$dist_c = d$Distance_m - mean(d$Distance_m)
    m = eDNA_lmer(Cq ~ dist_c + (1 + dist_c | FilterID), d,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5,
                  seed = 1, chains = 1, refresh = 0)

    rr = ranef(m)
    term = sub("\\[.*", "", rownames(rr))
    expect_equal(unique(term), m@random_groups)
    expect_false(is.unsorted(match(term, m@random_groups)))

    draws = unclass(m@fit$draws("rand_betas", format = "draws_matrix"))
    g = match(term, m@random_groups)
    # draws are stored to 6 significant digits, so allow relative rounding
    sums = rowsum(t(draws), g)
    expect_true(max(abs(sums) / rowsum(t(abs(draws)), g)) < 1e-5)

    res = summary(m)
    expect_true(all(paste0("sd(", m@random_groups, ")") %in% res$variable))
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

