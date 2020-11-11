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

pb_download("nets_mix_clean2k.RData",
            dest="data",
            tag="data.v.1")

load('data/nets_mix_clean2k.RData')

metric.net <- c('NODF',
                'H2',
                'vulnerability',
                'niche overlap')

#N = 1999
cores = 10


##Using the randomized networks built in 2_NetBuilding, calculate network level metrics
#NOTE: These steps are long, primarily due to the robustness simulations 
sexlvl <- mclapply(nets_mix_clean2k,
                   function(x) calcNets(x, metrics = metric.net),
                   mc.cores = cores)

save(sexlvl,file="data/sexlvl.RData")

pb_upload("data/sexlvl.RData",
          name="sexlvl.RData",
          tag = "data.v.1")

#robustness calculations. Species go extinct by degree, extinctions of female pollinators
#also result in the extinction of males of the same species. 
rob1a <- mclapply(nets_mix_clean2k,function(x){
  sites <- lapply(names(x), function(y){
    site.net <- x[[y]]
    extinct <- second.extinct1(site.net,nrep=100)
    site <- robustness1(extinct)
    site <- data.frame('lower'=site[[1]],'higher'=site[[2]],'SiteYr'=y)
    return(site)
  })
  sites <- do.call(rbind,sites)
  return(sites)
},mc.cores = cores)

save(rob1a,file='data/rob1a.RData')

pb_upload('data/rob1.RData',
            name='rob1a.RData',
            tag='data.v.1')



## ****************************************************************

pb_download("sexlvl.RData",
          dest="data",
          tag = "data.v.1")

pb_download('rob1.RData',
            dest='data',
            tag='data.v.1')

sexmet <- names(sexlvl1[[1]])[1:6]

#Compare observed metric values to simulated null network metric values
sexlvlProp50Z <- calcNullProp50(sexlvl,
                                sexmet,
                               zscore=T,
                               level="network")

robZ <- calcNullProp50(rob1a,
                       c('lower','higher'),
                       zscore=T,
                       level="network")

#add robustness zscores to the dataframe with the other metrics
sexlvlProp50Z$robustness.LL <- robZ$lower[match(sexlvlProp50Z$SpSiteYr,robZ$SpSiteYr)]
sexlvlProp50Z$robustness.HL <- robZ$higher[match(sexlvlProp50Z$SpSiteYr,robZ$SpSiteYr)]
sexlvlProp50Z <- sexlvlProp50Z[,c(1:6,8,9,7)]

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
sexmet <- names(sexlvl[[1]])[1:8]

overallTest(sexlvlProp50Z, sexmet, tails=2, zscore=T)

t.tester(sexlvlProp50Z,names(sexlvlProp50Z[,1:8]))
  #means for everything except for robustness are  significantly different from 0, 
  #all neg except for h2

