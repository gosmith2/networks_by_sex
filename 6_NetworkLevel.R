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

#THIS IS WHERE YOU TOGGLE INTERNAL NULLS. KEY: where sexlvlboot="y", and an epithet from a previous step (e.g., "boot5") = z
  #y = original: no bootstrapping, no internal nulls
  #yz = bootstrapped observations, no internal nulls
  #yzI = boostrapped observations, internal nulls, vaznull.fast
  #yzIB = boostrapped obs, internal nulls, basnull

sexlvlbootTIB.net <- mclapply(nets_mix_clean2kbootT,
                   function(x) calcNets(x, metrics = metric.net, null.sim=T,null.fun=basnull),
                   mc.cores = cores)


#sexlvlbootTI.f <- mclapply(c(1:2001),function(x){
#  df <- cbind(sexlvlbootTI[[x]],sexlvlbootTI.net[[x]][,c(1,2)])
#  df <- df[,c(6,7,1,2,3,4,5)]
#  return(df)
#},mc.cores=cores)

save(sexlvlbootTIB.net,file="data/sexlvlbootTIB.RData")

pb_upload("data/sexlvlboot.RData",
          name="sexlvlboot.RData",
          tag = "data.v.1")

#robustness calculations. Species go extinct by degree, extinctions of female pollinators
#also result in the extinction of males of the same species. 
rob1abootTI <- mclapply(nets_mix_clean2kbootT,function(x){
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

save(rob1abootST,file='data/rob1aboot5.RData')

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

sexmet <- names(sexlvlbootTIB[[1]])[1:6]

#Compare observed metric values to simulated null network metric values

sexlvlProp50ZbootTIB <- calcNullProp50(sexlvlbootTIB.net,
                                sexmet,
                               zscore=T,
                               level="network")

robZbootTI <- calcNullProp50(rob1abootTI,
                       c('lower','higher'),
                       zscore=T,
                       level="network")

#add robustness zscores to the dataframe with the other metrics
sexlvlProp50ZbootTIB$robustness.LL <- robZbootTI$lower[match(sexlvlProp50ZbootTIB$SpSiteYr,robZbootTI$SpSiteYr)]
sexlvlProp50ZbootTIB$robustness.HL <- robZbootTI$higher[match(sexlvlProp50ZbootTIB$SpSiteYr,robZbootTI$SpSiteYr)]
sexlvlProp50ZbootTIB <- sexlvlProp50ZbootTIB[,c(1:6,8,9,7)]

save(sexlvlProp50ZbootTIB,file="data/sexlvlProp50ZbootTIB.RData")
pb_upload('data/sexlvlProp50ZbootTIB.RData',
          name="sexlvlProp50ZbootTIB.RData",
          tag='data.v.1')

## ****************************************************************
pb_download("sexlvlProp50ZbootTIB.RData",
            dest="data",
            tag="data.v.1")


##Test: proportion of networks where values in the observed 
## network was significantly different than many of the simulations.
sexmet <- colnames(sexlvlProp50Zboot)[1:8]

overallTest(sexlvlProp50Z, sexmet, tails=2, zscore=T)

t.tester(sexlvlProp50Z,names(sexlvlProp50Z[,1:8]))
  #means for everything except for robustness are  significantly different from 0, 
  #all neg except for h2

