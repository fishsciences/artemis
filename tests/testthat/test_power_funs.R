context("Power functions")

vars = list(Intercept = -10.6,
            distance = c(0, 15, 50),
            volume = c(25, 50),
            biomass = 100,
            alive = 1,
            tech_rep = 1:10,
            rep = 1:3, Cq = 1)

X = expand.grid(vars)
betas = c(distance = -0.002, volume = 0.01, biomass = 1, alive = 1)

test_that("P-detect: basic tests", {
    expect_error(est_p_detect_lm(variable_levels = c(intecept = 1, Distance_m = 500),
                              betas = c(intercept = -10.6, Distance_m = -0.005, Volume_mL = 0.01),
                              ln_eDNA_sd = 1, 
                              n_rep = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5))
    
    expect_error(est_p_detect_lm(variable_levels = c(intecept = 1, Distance_m = 500, Volume_mL = 50),
                 betas = c(intercept = -10.6, Distance_m = -0.005),
                 ln_eDNA_sd = 1, 
                 n_rep = 6, std_curve_alpha = 21.2, std_curve_beta = -1.5))

    res = est_p_detect_lm(variable_levels = c(intecept = 1, Distance_m = 300, Volume_mL = 0),
                 betas = c(intercept = -10.6, Distance_m = -0.5, Volume_mL = -0.1),
                 ln_eDNA_sd = 1, 
                 n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    
    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)

    res = est_p_detect_lm(variable_levels = c(Distance_m = 50),
                 model_fit = ans, 
                 n_rep = 2)

    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == 1)

    res = est_p_detect_lm(variable_levels = c(Distance_m = 50),
                 model_fit = ans, 
                 n_rep = 2:6)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == length(2:6))


    expect_error(est_p_detect_lm(variable_levels = c(intecept = 1, Distance_m = 300,
                                             Volume_mL = c(0, 20)),
                              betas = c(intercept = -10.6, Distance_m = -0.5, Volume_mL = -0.1),
                              ln_eDNA_sd = 1, 
                              n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5))
   
    
})

test_that("P-detect: numeric", {
    # Mirrors analytic solution for estimating the detection probabilities
    # Here Volume is centered at 50ml and Biomass at 22.5
    known_betas = c(Intercept = -10.6, Distance_m = -0.009, Volume_mL = 0.008, N_fish = 0.023)
    std_curve = c(alpha = 21.073, beta = -1.453)
    upper_Cq = 40
    ln_conc_thresh = (upper_Cq - std_curve["alpha"]) / std_curve["beta"]

    max_dist = (ln_conc_thresh - known_betas["Intercept"]) / known_betas["Distance_m"]

    ans = est_p_detect_lm(variable_levels = c(Intercept = 1,
                                           Distance_m = max_dist,
                                           Volume_mL = 0, N_fish = 0),
                       betas = known_betas,
                       ln_eDNA_sd = 1, 
                       std_curve_alpha = std_curve["alpha"],
                       std_curve_beta = std_curve["beta"],
                       n_rep = 1)
    expect_true(all.equal(0.46, as.numeric(ans)))

    # Target of 5% detection
    tar_Cq = 42.46465 # 40 corresponds to the 5% quantile 
    tar_ln_conc = (tar_Cq - std_curve["alpha"]) / std_curve["beta"]

    max_dist = (tar_ln_conc - known_betas["Intercept"]) / known_betas["Distance_m"]

    ans = est_p_detect_lm(variable_levels = c(Intercept = 1,
                                           Distance_m = max_dist,
                                           Volume_mL = 0, N_fish = 0),
                       betas = known_betas,
                       ln_eDNA_sd = 1, 
                       std_curve_alpha = std_curve["alpha"],
                       std_curve_beta = std_curve["beta"],
                       n_rep = 1)

    #Should be close to 5%
    expect_true(abs(ans - 0.05) < 0.01)
})


test_that("Power functions", {
    skip_on_cran()
    v = list(Cq = 1, distance = c(0, 100, 300), rep = 1:3)
    fast = list(chains = 1, iter_warmup = 200, iter_sampling = 200,
                show_messages = FALSE, refresh = 0)

    ans = do.call(artemis:::est_power_lm,
                  c(list(Cq ~ distance, v, betas = c(-10.6, -0.01),
                         sigma_ln_eDNA = 1, std_curve_alpha = 21.2,
                         std_curve_beta = -1.5, n_sim = 2), fast))
    expect_equal(names(ans), c("(Intercept)", "distance"))
    expect_true(all(ans >= 0 & ans <= 1))

    ans = do.call(artemis:::est_power_lmer,
                  c(list(Cq ~ distance + (1|rep), v, betas = c(-10.6, -0.01),
                         sigma_ln_eDNA = 1, sigma_rand = 0.5,
                         std_curve_alpha = 21.2, std_curve_beta = -1.5,
                         type = "accuracy", n_sim = 2), fast))
    expect_equal(length(ans), 2L)

    ans = do.call(artemis:::est_power_range_lm,
                  c(list(Cq ~ distance, v, betas = c(-10.6, -0.01),
                         sigma_ln_eDNA = 1, std_curve_alpha = 21.2,
                         std_curve_beta = -1.5, rep_range = c(2, 4),
                         n_sim = 1), fast))
    expect_equal(names(ans), c("n_rep: 2", "n_rep: 4"))
})
