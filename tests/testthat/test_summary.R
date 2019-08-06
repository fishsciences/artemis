context("Simulating data")

vars = list(distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = .002, volume = -0.57, biomass = 1, alive = 1)

test_that("Summary methods: simulations", {
    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = 1, distance = 0.5, volume = 0),
                      sigma_Cq = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    res = summary(ans)

    expect_is(res, "eDNA_simulation.summary")
    expect_is(res, "data.frame")

    expect_true(nrow(res) == sum(sapply(ans@x, function(x) length(unique(x)))))
    expect_true(ncol(res) == 7)
    
    res2 = summary(ans, prob = 0.5)

    expect_true(ncol(res2) == 5)

    print(res)
})

test_that("Summary methods: model", {
    ans = eDNA_lm(Cq ~ Distance, eDNA_samples,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
    
    res = summary(ans)

    expect_is(res, "eDNA_model.summary")
    expect_is(res, "data.frame")

    expect_true(nrow(res) == ncol(ans@betas) + 1)
    expect_true(ncol(res) == 4)
    
    res2 = summary(ans, prob = 0.5)

    expect_true(ncol(res2) == 2)
    
    print(res)
    print(res2)
    
})