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
    ans = prep_data.sim(ml, 21, -1.5, 1, betas, 0.08, 
                    prior_int = normal(), prior_b = normal())

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_var_shared", "rand_sigma")], length) == 0))
    expect_true(all(sapply(ans[c("has_random", "n_rand")], `==`, 0)))

    ## lmer
    mler = gen_model_list_lmer(Cq ~ distance + (1|volume), X)
    ans = prep_data.sim(mler, 21, -1.5, 1, betas, 0.08, rand_sd = 0.1,
                    prior_int = normal(), prior_b = normal())

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_x", "rand_sigma")], length) != 0))
})
