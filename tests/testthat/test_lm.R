if(FALSE){
library(artemis)
library(elaphos)

# make sure model is up-to-date
compile_models("eDNA_lm.stan")

m = eDNA_lm(Cq ~ Distance_m, cvp02, 21, -1.5)
m
loo(m)

# only intercept
m1a = eDNA_lm(Cq ~ 1, cvp02, 21, -1.5)
m1a
loo(m1a)

# no intercept
m1b = eDNA_lm(Cq ~ -1 + Distance_m, cvp02, 21, -1.5)
m1b
loo(m1b)

# Default priors are not great - might look into autoscaling
m2a = eDNA_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5)
m2a
loo(m2a)

m2b = eDNA_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5,
              priors = normal(autoscale = FALSE))
m2b
loo(m2b)
# test that priors still work
m3 = eDNA_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5,
             priors = normal(1, c(0.02, 0.01)))
m3

compile_models("eDNA_lm_zinf.stan")
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5,
                            probability_zero = 0.08)
m4
loo(m4)

m5a = artemis:::eDNA_zinf_lm(Cq ~ 1, cvp02, 21, -1.5,
                            probability_zero = 0.08)
m5a
loo(m5a)
m5b = artemis:::eDNA_zinf_lm(Cq ~ -1 + Distance_m, cvp02, 21, -1.5,
                            probability_zero = 0.08)
m5b
loo(m5b)

}
