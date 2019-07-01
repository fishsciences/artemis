
                                        # Test for simulation functions
                                        # M. Espe
                                        # March 2019

context("Simulating data")
test_that("Simulate data", {

    vars = list(distance = c(0, 15, 50),
                volume = c(25, 50),
                biomass = 100,
                alive = 1,
                tech_rep = 1:10,
                rep = 1:3)

    X = expand.grid(vars)

    expect_equal(nrow(X), Reduce(`*`, sapply(vars, length)))

    betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

    rand = gen_rand_eff(X[,c("tech_rep", "rep")], c(.1, .1))
    expect_equal(nrow(rand), nrow(X))

    sim_data(vars, betas, c(0.1,0.1), .1)


    sim_data(X=X, betas = betas, sigma_rand =  c(0.1,0.1), sigma_Cq = .1)
})
