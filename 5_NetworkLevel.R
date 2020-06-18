#5: Network level analyses

##Loads dataframe of observed and simulated networks, calculates
## network-level metrics (including simulating extinction cascades), 
## and then compares observed to simulated values. 


library(piggyback)
library(fossil)
library(bipartite)
library(tidyverse) 
source('prepNets.R')

####----------------------------####
#### Prep for networklevel analyses
####----------------------------####

pb_download("nets_mix_clean.RData",
            dest="data",
            tag="data.v.1")

load('data/nets_mix_clean.RData')

metric.net <- c('NODF',
                'H2',
                'robustness',
                'vulnerability',
                'niche overlap')

N = 999
cores = 10


##Using the randomized networks built in 2_NetBuilding, calculate network level metrics
#NOTE: This is a long and computationally intensive step due to the robustness simulations 
sexlvl <- mclapply(nets_mix_clean,
                   function(x) calcNets(x, metrics = metric.net),
                   mc.cores = cores)
save(sexlvl,file="data/sexlvl.RData")

pb_upload("data/sexlvl.RData",
          name="sexlvl.RData",
          tag = "data.v.1")

## ****************************************************************

pb_download("sexlvl.RData",
          dest="data",
          tag = "data.v.1")

sexmet <- names(sexlvl[[1]])[1:8]

#Compare observed metric values to simulated null network metric values
sexlvlProp50 <- calcNullProp50(sexlvl,
                             sexmet,
                             zscore=F,
                             level="network")

#Repeat above to generate z-scores for plotting
sexlvlProp50Z <- calcNullProp50(sexlvl,
                                sexmet,
                               zscore=T,
                               level="network")

save(sexlvlProp50Z,file="data/sexlvlProp50Z.RData")
pb_upload('data/sexlvlProp50Z.RData',
          name="sexlvlProp50Z.RData",
          tag='data.v.1')

## ****************************************************************
pb_download("sexlvlProp50Z.RData",
            dest="data",
            tag="data.v.1")


##Test: proportion of networks where values in the observed 
## network was significantly different than many of the simulations.

overallTest(sexlvlProp50, sexmet, tails=2, zscore=F)

