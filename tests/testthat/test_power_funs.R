if(FALSE){
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

test_that("P-detect", {
    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 500),
                              betas = c(intercept = -10.6, Distance = -0.005, Volume = 0.01),
                              Cq_sd = 1, n_rep = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5))
    
    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 500, Volume = 50),
                 betas = c(intercept = -10.6, Distance = -0.005),
                 Cq_sd = 1, n_rep = 6, std_curve_alpha = 21.2, std_curve_beta = -1.5))

    res = est_p_detect(var_levels = c(intecept = 1, Distance = 300, Volume = 0),
                 betas = c(intercept = -10.6, Distance = -0.5, Volume = -0.1),
                 Cq_sd = 1, n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5)

    ans = eDNA_lm(Cq ~ Distance, eDNA_data,
                  std_curve_alpha = 21.2, std_curve_beta = -1.5)


    res = est_p_detect(var_levels = c(intecept = 1, Distance = 50),
                 betas = NULL, model_fit = ans, 
                 Cq_sd = 1, n_rep = 2, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == 1)

    res = est_p_detect(var_levels = c(intecept = 1, Distance = 50),
                 betas = NULL, model_fit = ans, 
                 Cq_sd = 1, n_rep = 2:6, std_curve_alpha = 21.2, std_curve_beta = -1.5)
    expect_true(nrow(res) == nrow(ans@betas))
    expect_true(ncol(res) == length(2:6))


    expect_error(est_p_detect(var_levels = c(intecept = 1, Distance = 300,
                                             Volume = c(0, 20)),
                              betas = c(intercept = -10.6, Distance = -0.5, Volume = -0.1),
                              Cq_sd = 1, n_rep = 3, std_curve_alpha = 21.2, std_curve_beta = -1.5))
   
    
})

test_that("Power function", {
    ans = est_power_lm(Cq ~ 1 + distance,
                    vars_list = list(Cq = 1, distance = c(0,100,300), rep = 1:3),
                    betas = c(intercept = -10.6, distance = -0.05),
                    sigma_Cq = 1, sigma_rep = 0.5,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 5,
                    type = "exclude_zero")

    ans = est_power_lm(Cq ~ 1 + distance,
                    vars_list = list(Cq = 1, distance = c(0,100,300), rep = 1:20),
                    betas = c(intercept = -10.6, distance = -0.05),
                    sigma_Cq = 1, sigma_rep = 0.5,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 5,
                    type = "exclude_zero")

    lapply(seq(2, 20, 2), function(i){
    ans = est_power_lm(Cq ~ 1 + distance,
                    vars_list = list(Cq = 1, distance = c(0,50,75,100,300), rep = 1:i),
                    betas = c(intercept = -10.6, distance = -0.04),
                    sigma_Cq = 1, sigma_rep = 0.5,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5, n_sim = 200,
                    type = "accuracy")
    })
 

})

}
