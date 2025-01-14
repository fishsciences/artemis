if(FALSE){
library(artemis)
library(tools)

compile_models("eDNA_lm_zinf.stan")
assertError(artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L))

m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | Volume_mL, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)

# should still work - drops the predictor to 0
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)

m4
loo(m4)

m5a = artemis:::eDNA_zinf_lm(Cq ~ 1, eDNA_data, 21, -1.5,
                             parallel_chains = 4L)
m5a
loo(m5a)
m5b = artemis:::eDNA_zinf_lm(Cq ~ -1 + Distance_m, eDNA_data, 21, -1.5,
                            probability_zero = 0.08)
m5b
loo(m5b)


x = formula(Cq ~ 1 + Distance_m | Distance_m + Volumn_mL, eDNA_data)


}
