context("Print methods")

## Figure out how to test this


vars = list(distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

test_that("Print methods: simulations", {
    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = 1, distance = 0.5, volume = 0),
                      sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    ## should not throw an error
    expect_is(print(ans), "eDNA_simulation")
})

test_that("Print methods: simulations", {
    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = 1, distance = 0.5, volume = 0),
                      sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    ## should not throw an error
    p_detect = est_p_detect(variable_levels = c(Intercept = 1, Distance = 100),
                        betas = c(Intercept = -10.5, Distance = -0.04),
                        ln_eDNA_sd = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5,
                        n_rep = 1:12)
    print(p_detect)

    expect_is(print(p_detect), "eDNA_p_detect")
})
