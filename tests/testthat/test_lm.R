library(artemis)
library(elaphos)

# make sure model is up-to-date
compile_models("eDNA_lm.stan")

m = eDNA_lm(Cq ~ Distance_m, cvp02, 21, -1.5)

m
