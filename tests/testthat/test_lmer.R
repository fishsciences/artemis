if(FALSE){
library(artemis)
## library(elaphos)

# make sure model is up-to-date
compile_models()

m = eDNA_lmer(Cq ~ Distance_m + (1|FilterID), eDNA_data, 21, -1.5)
m
loo(m)
ranef(m)

# slope
m0 = eDNA_lmer(Cq ~ Distance_m + (1 + Distance_m|FilterID), eDNA_data, 21, -1.5)
m0
loo(m)
ranef(m)

# only intercept
m1a = eDNA_lmer(Cq ~ 1 + (1|FilterID), eDNA_data, 21, -1.5)
m1a
loo(m1a)
ranef(m1a)

# no intercept

m1b = eDNA_lmer(Cq ~ -1 + Distance_m + (1|FilterID), eDNA_data, 21, -1.5)
m1b
loo(m1b)

# Default priors are not great - might look into autoscaling
m2a = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5)
m2a
loo(m2a)

m2b = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
              priors = normal(autoscale = FALSE))
m2b
loo(m2b)
# test that priors still work
m3 = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
             priors = normal(1, c(0.02, 0.01)))
m3

compile_models("eDNA_lm_zinf.stan")
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5)
m4
loo(m4)

m5a = artemis:::eDNA_zinf_lm(Cq ~ 1, eDNA_data, 21, -1.5)
m5a
loo(m5a)
m5b = artemis:::eDNA_zinf_lm(Cq ~ -1 + Distance_m, eDNA_data, 21, -1.5)
m5b
loo(m5b)

}

if(FALSE){
  library(artemis)
  ## Create count column - from Gblock std. curves
  # cp = 12.9 - 0.367 * Cq
  eDNA_data$count = 12.9 - 0.367 * eDNA_data$Cq
  eDNA_data$count = ceiling(eDNA_data$count)
  eDNA_data$count[eDNA_data$Cq == 40] = 0
  eDNA_data$count[eDNA_data$count < 0] = 1
  
  m = artemis:::eDNA_count_lm(count ~ 1, eDNA_data)

  m2 = artemis:::eDNA_count_lmer(count ~ 1 + (1|FilterID), eDNA_data)
  print(m2)

  ## Check estimates vs. rstanarm
  library(rstanarm)

  w = stan_glm(count ~ 1, eDNA_data, family = "poisson")

  ## similar estimates
  print(w, digits = 3)
  print(m)

  w2 = stan_glmer(count ~ 1 + (1|FilterID), eDNA_data, family = "poisson")

  ## Similar as well - need to get a more definitive test for this
  print(w2, digits = 3)
  print(m2)
}
