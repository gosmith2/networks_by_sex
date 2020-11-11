## 4: Dissimilarity indexing

## Loads networks and calculates the Morisita-Horn dissimilarity
## index for males and females within each species within each network.
## Lastly, the number of observations showing large differences from 
## the null is calculated. 

library(parallel)
library(piggyback)
library(vegan)
library(tidyverse) 
source('prepNets.R')


#Load clean networks
pb_download('nets_mix_clean2k.RData',
            dest="data",
            tag="data.v.1")

load("data/nets_mix_clean2kRdata")

#calculate the M-H index for species with at least 5 males and 5 females
distValues5.df <- distComp(nets_mix_clean2k,"horn",indiv=5)

#add siteyr column, remove NAs from sites that did not have species with the 
#minimum abundances
distValues5.df$SpSiteYr <- paste(distValues5.df$GenusSpecies,
                                 distValues5.df$SiteYr, sep="_")
distValues5.df <- filter(distValues5.df,distance!='NA')


#Save and upload distances
save(distValues5.df,file='data/distValues5.RData')

pb_upload("data/distValues5.RData",
          name="distValues5.RData",
          tag="data.v.1")

## ****************************************************************
#download to the previous skip step
pb_download("distValues5.RData",
            dest="data",
            tag="data.v.1")

#compare observed M-H distances to simulated null networks, generating zscores
diffDist5Zscore <- calcDistZ(distValues5.df,"SpSiteYr",zscore=T)
save(diffDist5Zscore,file='data/diffDist5Zscore.RData')

#diffDist5 <- calcDistZ(distValues5.df,"SpSiteYr",zscore=F)


## Test: proportion of species+sites where M-H distance in
## observed network was different than many of the simulations. The output
## is the proportion of observations (Sp+Site+Year) where the distance
## diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

#overallTest(diffDist5,"distanceZ",zscore=F,tails=2)
overallTest(diffDist5Zscore,"distanceZ",zscore=T,tails=2)

t.tester(diffDist5Zscore,'distanceZ')
  #mean for distribution of zscore values diff from 0 
  
#genera where a difference in M-H was observed
unique(word(unique(diffDist5Zscore$Level[abs(diffDist5Zscore$distanceZ)>1.96]),1))
  #10 genera
  
pb_upload("data/diffDist5Zscore.RData",
          name="diffDist5Zscore.RData",
          tag="data.v.1")



######-------Section S1: average together sites within a given dataset

load('data/zscore50_5.RData')
load('data/diffDist5Zscore.RData')

#network metrics
zscore50_5.avged <- zscore_avger(zscore50_5.df,metric.ls)

overallTest(zscore50_5.avged, metric.ls, tails =2,zscore=T)
t.tester(zscore50_5.avged,metric.ls)

#M-H distance
diffDist5Zscore.avged <- distZ_avger(diffDist5Zscore)

overallTest(diffDist5Zscore.avged,"distanceZ",zscore=T,tails=2)
t.tester(diffDist5Zscore.avged,'distanceZ')

