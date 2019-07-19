
                                        # Test for simulation functions
                                        # M. Espe
                                        # March 2019

context("Simulating data")

vars = list(distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

test_that("Simulate data: lm", {

    expect_equal(nrow(X), Reduce(`*`, sapply(vars, length)))

    betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = 1, distance = 0.5)))

    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = 1, distance = 0.5)))

    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = 1, distance = 0.5, volume = 0),
                      sigma_Cq = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(ans, "eDNA_simulation")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
    
})

test_that("Simulate data: lmer", {

    #rand effects not provided
    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                               vars))

    ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                        vars,
                        betas = c(intercept = 1, distance = 0.5, volume = 0),
                        sigma_Cq = 1,
                        sigma_rand = 0.1,
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5)
    expect_is(ans, "eDNA_simulation")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
})

test_that("Prep data", {
    ## lm
    ml = gen_model_list_lm(Cq ~ distance, X)
    ans = prep_sim(ml, 21, -1.5, 1, betas)

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_var_shared", "rand_sigma")], length) == 0))
    expect_true(all(sapply(ans[c("has_rand", "n_rand_var", "n_rand_total")], `==`, 0)))

    ## lmer
    mler = gen_model_list_lmer(Cq ~ distance + (1|volume), X)
    ans = prep_sim(mler, 21, -1.5, 1, betas, 0.1)

    expect_is(ans, "list")

    expect_true(all(sapply(ans[c("groups", "rand_var_shared", "rand_sigma")], length) != 0))
    expect_true(all(sapply(ans[c("has_rand", "n_rand_var", "n_rand_total")], `!=`, 0)))
})

test_that("Gen model list", {
    ml = gen_model_list_lm(Cq ~ distance, X)
    
    expect_is(ml, "list")
    expect_true(all(names(ml) %in% c("x","y")))

    mler = gen_model_list_lmer(Cq ~ distance + (1|volume), X)
    expect_is(mler, "list")
    expect_true(all(names(ml) %in% c("x","y", "groups", "n_grp", "n_levels")))
    expect_true(ncol(mler$groups) == mler$n_grp)
    expect_true(length(mler$n_levels) == ncol(mler$groups))
    expect_true(all(sapply(mler$groups, max) == mler$n_levels))
   
})

test_that("Multiple groups", {
    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                        vars,
                        betas = c(intercept = 1, distance = 0.5, volume = 0),
                        sigma_Cq = 1,
                        sigma_rand = c(0.1),
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5))

    ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                        vars,
                        betas = c(intercept = 1, distance = 0.5, volume = 0),
                        sigma_Cq = 1,
                        sigma_rand = c(0.1, 0.1),
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5)


})
