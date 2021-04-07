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
    expect_error(est_p_detect(variable_levels = c(intecept = 1, Distance_m = 500),
                              betas = c(intercept = -10.6, Distance_m = -0.005, Volume_mL = 0.01),
                              ln_eDNA_sd = 1, 
                              n_rep = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5))
    
    expect_error(est_p_detect(variable_levels = c(intecept = 1, Distance_m = 500, Volume_mL = 50),
                 betas = c(intercept = -10.6, Distance_m = -0.005),
                 ln_eDNA_sd = 1, 
                 n_rep = 6, std_curve_alpha = 21.2, std_curve_beta = -1.5))

    res = est_p_detect(variable_levels = c(intecept = 1, Distance_m = 300, Volume_mL = 0),
                 betas = c(intercept = -10.6, Distance_m = -0.5, Volume_mL = -0.1),
                 ln_eDNA_sd = 1, 
                 n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    
    ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)

    res = est_p_detect(variable_levels = c(Distance_m = 50),
                 model_fit = ans, 
                 n_rep = 2)

    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == 1)

    res = est_p_detect(variable_levels = c(Distance_m = 50),
                 model_fit = ans, 
                 n_rep = 2:6)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == length(2:6))


    expect_error(est_p_detect(variable_levels = c(intecept = 1, Distance_m = 300,
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

    ans = est_p_detect(variable_levels = c(Intercept = 1,
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

    ans = est_p_detect(variable_levels = c(Intercept = 1,
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


if(FALSE){ #These can take a really long time to run
    
test_that("Power function", {
    ans = est_power_lm(Cq ~ 1 + distance,
                       variable_list = list(Cq = 1, distance = c(0,100,300), rep = 1:3),
                       betas = c(intercept = -10.6, distance = -0.05),
                       sigma_Cq = 1, sigma_rep = 0.5,
                       std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 5,
                       type = "exclude_zero")

    ans = est_power_lm(Cq ~ 1 + distance,
                    variable_list = list(Cq = 1, distance = c(0,100,300), rep = 1:20),
                    betas = c(intercept = -10.6, distance = -0.05),
                    sigma_Cq = 1, sigma_rep = 0.5,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 5,
                    type = "exclude_zero")

    lapply(seq(2, 20, 2), function(i){
    ans = est_power_lm(Cq ~ 1 + distance,
                    variable_list = list(Cq = 1, distance = c(0,50,75,100,300), rep = 1:i),
                    betas = c(intercept = -10.6, distance = -0.04),
                    sigma_Cq = 1, sigma_rep = 0.5,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 200,
                    type = "accuracy")
    })
 

})

}
