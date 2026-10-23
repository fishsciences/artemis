context("Prep data")

vars = list(distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

test_that("Prep data", {
    ## lm
    ml = gen_model_list_lm(Cq ~ distance, X)
    ans = prep_data.sim(ml, 21, -1.5, 1, betas, 0.08)

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_var_shared", "rand_sigma")], length) == 0))
    expect_true(all(sapply(ans[c("has_random", "n_rand")], `==`, 0)))

    ## lmer
    mler = gen_model_list_lmer(Cq ~ distance + (1|volume), X)
    ans = prep_data.sim(mler, 21, -1.5, 1, betas, 0.08, rand_sd = 0.1)

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_x", "rand_sigma")], length) != 0))
})

test_that("Count response conversion", {
    d = data.frame(y = c(0, 3, 10), x = 1:3)
    ml = gen_model_list_lm(y ~ x, d)
    # integer-valued doubles convert silently
    ans = expect_silent(prep_data.count(ml, prior_int = normal(), prior_b = normal()))
    expect_true(is.integer(ans$y))

    ml$y = c(0.4, 3, 10.6)
    expect_warning(prep_data.count(ml, prior_int = normal(), prior_b = normal()),
                   "not integer")
    ans = suppressWarnings(prep_data.count(ml, prior_int = normal(), prior_b = normal()))
    expect_equal(ans$y, c(0L, 3L, 11L))
})

test_that("Random effect matrices stay matrices with one censored row", {
    d = eDNA_data[1:30, ]
    d$Cq = pmin(d$Cq, 39)
    d$Cq[1] = 40
    ml = gen_model_list_lmer(Cq ~ Distance_m + (1|FilterID), d)
    md = prep_data.lmer(ml, 21.2, -1.5, 0, 40, normal(), normal(), exponential())
    expect_true(is.matrix(md$X_cens_r))
    expect_equal(dim(md$X_cens_r), c(1L, md$K_r))
})

test_that("remove_intercept only drops the intercept column", {
    x = cbind("(Intercept)" = 1, Intercept_dist = 1:3, b = 3:1)
    expect_equal(colnames(remove_intercept(x)), c("Intercept_dist", "b"))
})

test_that("Autoscaled priors shrink for wide predictors", {
    x = data.frame(a = c(0, 100, 200, 300), b = c(0, 1, 0, 1))
    y = c(1, 2, 3, 5)
    p = prep_priors(normal(0, 1), x, y)
    expect_equal(as.vector(p$scale), c(sd(y) / sd(x$a), sd(y)))
})
