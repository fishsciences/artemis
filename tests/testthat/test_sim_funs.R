
                                        # Test for simulation functions
                                        # M. Espe
                                        # March 2019

context("Simulating data")

vars = list(Intercept = -10.6,
            distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = -0.002, volume = 0.01, biomass = 1, alive = 1)

test_that("Simulate data: lm", {
    expect_error(sim_eDNA_lm( ~ distance + volume, vars,
                             betas = c(intercept = -10.6, distance = -0.05)))
    
    expect_equal(nrow(X), Reduce(`*`, sapply(vars, length)))

    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = -10.6, distance = -0.05)))
    
    expect_error(sim_eDNA_lm(Cq ~ distance + volume, vars,
                             betas = c(intercept = -10.6, distance = -0.05)))

    ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                      betas = c(intercept = -10.6, distance = -0.05, volume = 0.1),
                      sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_is(ans, "eDNA_simulation")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
    
})

test_that("Simulate data: lmer", {
    expect_error(sim_eDNA_lmer( ~ distance + volume + (1|rep),
                               vars))

    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                               vars,
                               betas = c(intercept = -10.6, distance = -0.05,
                                         volume = 0.01),
                               sigma_ln_eDNA = 1,
                               std_curve_alpha = 21.2,
                               std_curve_beta = -1.5))

    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                               vars,
                               betas = c(intercept = -10.6, volume = 0.01),
                               sigma_ln_eDNA = 1,
                               sigma_rand = 0.1,
                               std_curve_alpha = 21.2,
                               std_curve_beta = -1.5))

    #rand effects not provided
    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                               vars))

    ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep),
                        vars,
                        betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),
                        sigma_ln_eDNA = 1,
                        sigma_rand = 0.1,
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5)
    expect_is(ans, "eDNA_simulation")
    ## expect_true(all(names(ans) %in% c("x", "Cq_star", "ln_conc")))
})


test_that("Gen model list", {
    ml = gen_model_list_lm(Cq ~ distance, X)
    
    expect_is(ml, "list")
    expect_true(all(names(ml) %in% c("x","y")))

    mler = gen_model_list_lmer(Cq ~ distance + (1|volume), X)
    expect_is(mler, "list")
    expect_true(all(names(ml) %in% c("x","y", "groups", "n_grp", "n_levels")))
    expect_true(ncol(mler$rand_x) == mler$n_rand)
    expect_true(all(sapply(mler$groups, max) == mler$n_levels))
   
})

test_that("Multiple groups", {
    expect_error(sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                        vars,
                        betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),
                        sigma_ln_eDNA = 1,
                        sigma_rand = c(0.1),
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5))

    ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                        vars,
                        betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),
                        sigma_ln_eDNA = 1,
                        sigma_rand = c(0.1, 0.1),
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5)

    

})

test_that("Multiple sims", {

    ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                        vars, n_sim = 10,
                        betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),
                        sigma_ln_eDNA = 1,
                        sigma_rand = c(0.1, 0.1),
                        std_curve_alpha = 21.2,
                        std_curve_beta = -1.5)
    expect_is(ans, "eDNA_simulation")
    

})

# might want to move these eventually to a methods test file
test_that("Methods", {
  ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
                    betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),,
                    sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)


  print(ans)
  expect_is(print(ans), "eDNA_simulation")

  a = summary(ans)

  expect_is(a, "eDNA_simulation.summary")
  expect_is(a, "data.frame")
  
  ans2 = as(ans, "data.frame")

  expect_is(ans2, "data.frame")
  ans3 = as.data.frame(ans)

  expect_is(ans3, "data.frame")
  expect_equal(ans2, ans3)

  ans = sim_eDNA_lmer(Cq ~ distance + volume + (1|rep) + (1|tech_rep),
                      vars, n_sim = 10,
                      betas = c(intercept = -10.6, distance = -0.05, volume = 0.01),
                      sigma_ln_eDNA = 1,
                      sigma_rand = c(0.1, 0.1),
                      std_curve_alpha = 21.2,
                      std_curve_beta = -1.5)

  print(ans)
  a = summary(ans)

  expect_is(a, "eDNA_simulation.summary")
  expect_is(a, "data.frame")

})

test_that("Simulation accuracy", {
    ans = sim_eDNA_lm(Cq ~ 1, vars,
                      betas = c(intercept = -15),
                      sigma_ln_eDNA = 1e-5, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    
    expect_is(ans, "eDNA_simulation")
    expect_true(all(round(ans@ln_conc_hat) == -15))
    expect_true(all(ans@Cq_star == 40))

    ans = sim_eDNA_lm(Cq ~ 1, vars,
                      betas = c(intercept = -10),
                      sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    expect_true(abs(sd(ans@Cq_star) - 1) < 1.0)
    
    ans = sim_eDNA_lm(Cq ~ 1 + distance, vars,
                      betas = c(intercept = -10.6, distance = -11),
                      sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    tt = tapply(ans@Cq_star, ans@x$distance, unique)
    expect_true(all(tt[c("15", "50") == 40]))
    
    tt = tapply(ans@Cq_star, ans@x$distance, sd)
    expect_true(all(abs(tt[1] - 1) < 1)) # the sd of Cq should be close to 1
})
