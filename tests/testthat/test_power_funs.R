context("Power calculations")
test_that("Test Power", {

    vars = list(intercept = 1,
                distance = c(0, 100, 500),
                volume = c(-10, 0, 10),
                tech_rep = 1:10,
                rep = 1:3)
            
    betas = c(intercept = -10, distance = .02, volume = -0.057) 

    model_data = prep_model_data(vars, betas, sigma_rand = c(0.1, 0.1),
                                 sigma_Cq = 0.1)

    mod = stanmodels$eDNA_power
    
    fit = sampling(mod, model_data, iter = 400, cores = 1)
    
# Default test
ans = est_power(mod, 5, vars_list = vars, betas = betas,
                sigma_rand = c(0.1, 0.1), sigma_Cq = 0.1)
ans

# User supplied matrix
X = expand.grid(vars)

X$bob = X$distance * X$volume

ans = est_power(mod, 1, mod_matrix = X, betas = c(betas, bob = .1),
                sigma_rand = c(0.1, 0.1), sigma_Cq = 0.1)

# Interval test rather than power
ans = est_power(mod, 5, vars_list = vars, betas = betas,
                sigma_rand = c(0.1, 0.1), sigma_Cq = 0.1, result_type = "interval")

ans
})
