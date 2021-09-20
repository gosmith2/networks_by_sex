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

pb_download("nets_mix_clean10kO.RData",
            dest="data",
            tag="data.v.1")
load('data/nets_mix_clean10kO.RData')

pb_download("nets_mix_clean10kS.RData",
            dest="data",
            tag="data.v.1")

load('data/nets_mix_clean10kS.RData')


metric.net <- c('NODF',
                'H2',
                'vulnerability')

cores = 10


##Using the randomized networks built in 2_NetBuilding, calculate network level metrics

sexlvlO.net <- mclapply(nets_mix_clean10kO,
                   function(x) calcNets(x, metrics = metric.net),
                   mc.cores = cores)
sexlvlS.net <- mclapply(nets_mix_clean10kS,
                        function(x) calcNets(x, metrics = metric.net),
                        mc.cores = cores)


save(sexlvlS.net,file="data/sexlvlS.RData")
save(sexlvlO.net,file="data/sexlvlO.RData")


#robustness calculations. Species go extinct by degree, extinctions of female pollinators
#also result in the extinction of males of the same species. 
rob1O <- mclapply(nets_mix_clean10kO,function(x){
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

rob1S <- mclapply(nets_mix_clean10kS,function(x){
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

save(rob1O,file='data/rob1O.RData')
save(rob1S,file='data/rob1S.RData')




## ****************************************************************

sexmet <- names(sexlvlO.net[[1]])[1:4]

#Compare observed metric values to simulated null network metric values

sexlvlProp50ZO <- calcNullProp50(sexlvlO.net,
                                sexmet,
                               zscore=T,
                               level="network")
sexlvlProp50ZS <- calcNullProp50(sexlvlS.net,
                                 sexmet,
                                 zscore=T,
                                 level="network")

robZO <- calcNullProp50(rob1O,
                       c('lower','higher'),
                       zscore=T,
                       level="network")
robZS <- calcNullProp50(rob1S,
                        c('lower','higher'),
                        zscore=T,
                        level="network")

#add robustness zscores to the dataframe with the other metrics
sexlvlProp50ZO$robustness.LL <- robZO$lower[match(sexlvlProp50ZO$SpSiteYr,robZO$SpSiteYr)]
sexlvlProp50ZO$robustness.HL <- robZO$higher[match(sexlvlProp50ZO$SpSiteYr,robZO$SpSiteYr)]
sexlvlProp50ZO <- sexlvlProp50ZO[,c(1:4,6,7,5)]

sexlvlProp50ZS$robustness.LL <- robZS$lower[match(sexlvlProp50ZS$SpSiteYr,robZS$SpSiteYr)]
sexlvlProp50ZS$robustness.HL <- robZS$higher[match(sexlvlProp50ZS$SpSiteYr,robZS$SpSiteYr)]
sexlvlProp50ZS <- sexlvlProp50ZS[,c(1:4,6,7,5)]

save(sexlvlProp50ZO,file="data/sexlvlProp50ZO.RData")
save(sexlvlProp50ZS,file="data/sexlvlProp50ZS.RData")

#pb_upload('data/sexlvlProp50ZS.RData',
#          name="sexlvlProp50ZS.RData",
#          tag='data.v.1')
#pb_upload('data/sexlvlProp50ZO.RData',
#          name="sexlvlProp50ZO.RData",
#          tag='data.v.1')


## ****************************************************************
#pb_download("sexlvlProp50ZS.RData",
#            dest="data",
#            tag="data.v.1")
#pb_download("sexlvlProp50ZO.RData",
#            dest="data",
#            tag="data.v.1")
load('data/sexlvlProp50ZS.RData')
load('data/sexlvlProp50ZO.RData')

#add columns with other variables
sexmet <- colnames(sexlvlProp50ZO)[1:6]
zNetO <- zscore_net_joiner(sexlvlProp50ZO,nets_mix_clean10kO,spec_all,"O")
zNetS <- zscore_net_joiner(sexlvlProp50ZS,nets_mix_clean10kS,spec_all,"S")

save(zNetO,file='data/zNetO.RData')
save(zNetS,file='data/zNetS.RData')

pb_upload('data/zNetO.RData',
          name="zNetO.RData",
          tag='data.v.1')
pb_upload('data/zNetS.RData',
          name="zNetS.RData",
          tag='data.v.1')

##Test: Do network-level metrics differ between the observed and simulated networks?

nodf.O <- lm(NODF~dataset,data=zNetO)
summary(nodf.O)
h2.O <- lm(H2~dataset,data=zNetO)
summary(h2.O)
gen.O <- lm(generality.HL~dataset,data=zNetO)
summary(gen.O)
vuln.O <- lm(vulnerability.LL~dataset,data=zNetO)
summary(vuln.O)
robH.O <- lm(robustness.HL~dataset,data=zNetO)
summary(robH.O)
robL.O <- lm(robustness.LL~dataset,data=zNetO)
summary(robL.O)



#--- Validate with "S" data
nodf.S <- lm(NODF~dataset,data=zNetS)
summary(nodf.S)
h2.S <- lm(H2~dataset,data=zNetS)
summary(h2.S)
gen.S <- lm(generality.HL~dataset,data=zNetS)
summary(gen.S)
vuln.S <- lm(vulnerability.LL~dataset,data=zNetS)
summary(vuln.S)
robH.S <- lm(robustness.HL~dataset,data=zNetS)
summary(robH.S)
robL.S <- lm(robustness.LL~dataset,data=zNetS)
summary(robL.S)

