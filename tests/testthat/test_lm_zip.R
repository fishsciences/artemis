if(FALSE){
library(artemis)
library(tools)

compile_models("eDNA_lm_zinf.stan")

##no zip - should error
assertError(artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L))

m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | Distance_m, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)
summary(plogis(m4@fit[,"nz_alpha[1]"]))

# should still work - drops the predict
m4 = artemis:::eDNA_zinf_lm(Cq ~ Distance_m + Volume_mL | 1, eDNA_data, 21, -1.5,
                                   parallel_chains = 4L)


## Need simulated dataset to test accuracy
df = data.frame(Cq = rnorm(200, 35, 1))
df$Cq[sample(1:200, 50)] = 40

m = artemis:::eDNA_zinf_lm(Cq ~ 1 | 1,df , 21, -1.5)
summary(plogis(m@fit[,"nz_alpha[1]"]))

}
