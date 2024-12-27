if(FALSE){
library(artemis)

# make sure model is up-to-date
compile_models("eDNA_lm.stan")

m = eDNA_lm(Cq ~ Distance_m, eDNA_data, 21, -1.5, parallel_chains = 4L)
m
loo(m)

# only intercept
m1a = eDNA_lm(Cq ~ 1, eDNA_data, 21, -1.5, parallel_chains = 4L)
m1a
loo(m1a)

# no intercept
m1b = eDNA_lm(Cq ~ -1 + Distance_m, eDNA_data, 21, -1.5, parallel_chains = 4L)
m1b
loo(m1b)

# Default priors are not great - might look into autoscaling
m2a = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5, parallel_chains = 4L)
m2a
loo(m2a)

m2b = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
              priors = normal(autoscale = FALSE), parallel_chains = 4L)
m2b
loo(m2b)
# test that priors still work
m3 = eDNA_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
             priors = normal(1, c(0.02, 0.01)), parallel_chains = 4L)
m3

## compile_models("eDNA_lm_zinf.stan")

m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5, parallel_chains = 4L)
m4
loo(m4)

m5a = artemis:::eDNA_zinf_lm(Cq ~ 1, eDNA_data, 21, -1.5, parallel_chains = 4L)
m5a
loo(m5a)
m5b = artemis:::eDNA_zinf_lm(Cq ~ -1 + Distance_m, eDNA_data, 21, -1.5, parallel_chains = 4L)
m5b
loo(m5b)

}
