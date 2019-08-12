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
                      sigma_Cq = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    ## should not throw an error
    expect_is(print(ans), "eDNA_simulation")
})
