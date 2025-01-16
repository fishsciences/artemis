if(FALSE){
library(artemis)
library(tools)

compile_models("eDNA_lm_zinf.stan")
assertError(artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L))

m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | Distance_m, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L, adapt_delta = 0.99 )

# should still work - drops the predict
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)


## Need simulated dataset to test accuracy
df = data.frame(Cq = rnorm(200, 35, 1))
df$Cq[sample(1:200, 50)] = 40

m = artemis:::eDNA_zinf_lm(Cq ~ 1 | 1,df , 21, -1.5)
summary(plogis(m@fit[,"z_alpha"]))

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
