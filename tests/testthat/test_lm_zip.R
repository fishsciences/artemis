if(FALSE){
library(artemis)
library(tools)

compile_models("eDNA_lm_zinf.stan")

##no zip - should error
assertError(eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L))

m4 = eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | Distance_m, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)
summary(plogis(m4@fit[,"nz_alpha[1]"]))

# should still work - drops the predict
m4 = eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)


## Need simulated dataset to test correctness
n = 1e5 ## need sufficiently large sample to get accurate estimates
df = data.frame(Cq = rnorm(n, 35, 1))
df$z = sample(c(TRUE, FALSE), nrow(df), replace = TRUE, prob = c(0.75, 0.25))
df$Cq[!df$z] = 40

## takes a long time due to data dims
m = eDNA_zinf_lm(Cq ~ 1 | 1,df , 21, -1.5, parallel_chains = 4L, adapt_delta = 0.99)
summary(plogis(m@fit$draws(variables = "nz_alpha[1]", format = "draws_df")[[1]]))


m = eDNA_zinf_lm(Cq ~ 1 | z,df , 21, -1.5, parallel_chains = 4L)
summary(plogis(m@fit$draws(variables = "nz_beta[1]", format = "draws_df")[[1]]))


compile_models("eDNA_lmer_zinf.stan")

m5 = eDNA_zinf_lmer(Cq ~ Distance_m + Volume_mL + (1|FilterID) | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L, adapt_delta = 0.99)

m6 = eDNA_zinf_lmer(Cq ~ Distance_m + (1|FilterID) | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L, adapt_delta = 0.99)
m6

}
