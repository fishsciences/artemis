
                                        # Test for simulation functions
                                        # M. Espe
                                        # March 2019

context("Simulating data")
test_that("Simulate data: lm", {

    vars = list(distance = c(0, 15, 50),
                volume = c(25, 50),
                biomass = 100,
                alive = 1,
                tech_rep = 1:10,
                rep = 1:3, Cq = 1)

    X = expand.grid(vars)

    expect_equal(nrow(X), Reduce(`*`, sapply(vars, length)))

    betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = 1, distance = 0.5)))

    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = 1, distance = 0.5)))

    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = 1, distance = 0.5, volume = 0),
                      sigma_Cq = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(ans, "list")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
    
})

test_that("Simulate data: lmer", {

    vars = list(distance = c(0, 15, 50),
                volume = c(25, 50),
                biomass = 100,
                alive = 1,
                tech_rep = 1:10,
                rep = 1:3, Cq = 1)

    X = expand.grid(vars)

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
    expect_is(ans, "list")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
})
