#--------------------------------------------#
# Setup Script for eDNA sample data analysis, 
# Data source: Delta (Miner Slough) Experiments 2-4
# Myfanwy Johnston, Feb 2021
# updated Thu Apr 15 12:41:09 2021 ------------------------------
#--------------------------------------------#
#--------------------------------------------#
if(FALSE){
d = rbind(elaphos::cvp01, elaphos::cvp02)
# Need standard curve
i = match(d$StdCrvID, elaphos::StdCrvKey$StdCrvID)
d$StdCrvAlpha_lnForm = elaphos::StdCrvKey$StdCrvAlpha_lnForm[i]
d$StdCrvBeta_lnForm = elaphos::StdCrvKey$StdCrvBeta_lnForm[i]

eDNA_data = d[ , c("Date", "FilterID", "TechRep", "Cq", "Distance_m", "Volume_mL", "Biomass_N", "StdCrvAlpha_lnForm", "StdCrvBeta_lnForm" )]

save(eDNA_data, file = "data/sample_data.rda")

#-------------------------------------------------------#










#-------------------------------------------------------#
# OLD artemis versions sample data - deprecated
#-------------------------------------------------------#

if(FALSE){library(dplyr)
# Experimental data

unique((read.csv("~/NonDropboxRepos/DS_predict/data/Exps2-4_20190124.csv", stringsAsFactors = FALSE) %>% 
                               select(-FlowMeterID, -Biomass))$Date)

cqdat = rbind(as.data.frame(read.csv("~/NonDropboxRepos/DS_predict/data/Exps2-4_20190124.csv", stringsAsFactors = FALSE) %>% 
                               select(-FlowMeterID, -Biomass) %>% 
                               mutate(Date = lubridate::ymd(Date))),
               as.data.frame(readxl::read_excel("~/NonDropboxRepos/DS_predict/data/Exp4.5_20190321.xlsx", sheet = 1) %>% 
                               select(-Biomass) %>% 
                               mutate(Date = lubridate::ymd(Date))
                             )
               )

str(cqdat)

cqdat = cqdat %>% 
  filter(Experiment != "Control" & Experiment != "Max volume" &
           Distance != "Control" & Distance != "Post") %>%  # filter out control samples
  filter(Volume < 180) %>% 
  mutate(Experiment = as.numeric(Experiment),
         Distance = as.numeric(Distance),
         Volume = as.numeric(Volume),
         SampleNumber = as.numeric(SampleNumber))

# Make sample data for package:
# Fri Jul  5 12:35:19 2019 ------------------------------

cqdat %>%
  filter(Volume !=80) -> smd


table(smd$Experiment, smd$SampleNumber)
 
table(smd$Distance, smd$Volume)

smd %>%
  arrange(Experiment, Volume, Distance, SampleNumber) %>%
  mutate(SampleID = data.table::rleid(SampleNumber)) -> smd

table(smd$SampleID)
table(smd$Experiment)
smd <- select(smd, -SampleTime, -Experiment, FilterNumber = SampleNumber)
 
ggplot(smd, aes(x = Volume, y = Cq)) +
  geom_point() +
  facet_grid(~Distance)

unique(smd$Date)
# unique(cqdat$Date)
 
# eDNA_samples <- select(smd, Date, SampleID, TechnicalRep, FilterNumber, Distance, Volume, Cq)
# #save(eDNA_samples, file = "../data/sample_data.rda")
# ggplot(eDNA_samples) +
#   geom_jitter(aes(x = Distance, y = Cq, color = factor(Volume)))

}
}
