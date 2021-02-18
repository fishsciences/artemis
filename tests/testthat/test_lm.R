library(artemis)
library(elaphos)

# make sure model is up-to-date
compile_models("eDNA_lm.stan")

m = eDNA_lm(Cq ~ Distance_m, cvp02, 21, -1.5)

m

m2 = eDNA_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5)
m2

# test that priors still work
m3 = eDNA_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5,
             priors = normal(1, c(0.02, 0.01)))
m3

compile_models("eDNA_lm_zinf.stan")
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, cvp02, 21, -1.5,
                            probability_zero = 0.08)
m4
