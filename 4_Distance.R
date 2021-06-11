## 4: Dissimilarity indexing

## Loads networks and calculates the Morisita-Horn dissimilarity
## index for males and females within each species within each network.
## Lastly, the number of observations showing large differences from 
## the null is calculated. 

library(parallel)
library(piggyback)
library(vegan)
library(tidyverse) 
library(stringr)
source('prepNets.R')


#Load clean networks
pb_download('nets_mix_clean2k.RData',
            dest="data",
            tag="data.v.1")

load("data/nets_mix_clean1kboot.Rdata")

#calculate the M-H index for species with at least 5 males and 5 females
distValues5S.df <- distComp(nets_mix_clean10kS,"horn",indiv=5)
distValues5O.df <- distComp(nets_mix_clean10kO,"horn",indiv=5)


#add siteyr column, remove NAs from sites that did not have species with the 
#minimum abundances
distValues5O.df$SpSiteYr <- paste(distValues5O.df$GenusSpecies,
                                 gsub("\\.","_",distValues5O.df$SiteYr), sep="_")
distValues5O.df <- filter(distValues5O.df,distance!='NA')


distValues5S.df$SpSiteYr <- paste(distValues5S.df$GenusSpecies,
                                       gsub("\\.","_",distValues5S.df$SiteYr), sep="_")
distValues5S.df <- filter(distValues5S.df,distance!='NA')



#Save and upload distances
save(distValues5O.df,file='data/distValues5O.RData')
save(distValues5S.df,file='data/distValues5S.RData')


pb_upload("data/distValues5S.RData",
          name="distValues5S.RData",
          tag="data.v.1")

## ****************************************************************
#download to the previous skip step
pb_download("distValues5O.RData",
            dest="data",
            tag="data.v.1")

#compare observed M-H distances to simulated null networks, generating zscores
diffDist5ZscoreO <- calcDistZ(distValues5O.df,"SpSiteYr",zscore=T)
diffDist5ZscoreS <- calcDistZ(distValues5S.df,"SpSiteYr",zscore=T)

save(diffDist5ZscoreS,file='data/diffDist5ZscoreS.RData')
save(diffDist5ZscoreO,file='data/diffDist5ZscoreO.RData')

pb_upload("data/diffDist5ZscoreS.RData",
          name="diffDist5ZscoreS.RData",
          tag="data.v.1")
pb_upload("data/diffDist5ZscoreO.RData",
          name="diffDist5ZscoreO.RData",
          tag="data.v.1")

#diffDist5 <- calcDistZ(distValues5.df,"SpSiteYr",zscore=F)


## Test: proportion of species+sites where M-H distance in
## observed network was different than many of the simulations. The output
## is the proportion of observations (Sp+Site+Year) where the distance
## diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

#overallTest(diffDist5,"distanceZ",zscore=F,tails=2)
overallTest(diffDist5ZscoreO,"distanceZ",zscore=T,tails=2)

t.tester(diffDist5ZscoreO,'distanceZ')
  #mean for distribution of zscore values diff from 0 
  
#genera where a difference in M-H was observed
unique(word(unique(diffDist5Zscoreboot$Level[abs(diffDist5Zscore$distanceZ)>1.96]),1))
  #10 genera
  #with ST, there are 27 genera
  
pb_upload("data/diffDist5ZscorebootST.RData",
          name="diffDist5ZscorebootST.RData",
          tag="data.v.1")

pb_download('diffDist5ZscoreO.RData',dest='data',tag='data.v.1')

######-------Section S1: average together sites within a given dataset

load('data/zscore50_5.RData')
load('data/diffDist5ZscoreS.RData')
load('data/spec_all.RData')
metric.ls <- c("degree","weighted.betweenness",
               "weighted.closeness","d")

#network metrics
zscore50_5.avged <- zscore_avger(zscore50_5.df,metric.ls)

overallTest(zscore50_5.avged, metric.ls, tails =2,zscore=T)
t.tester(zscore50_5.avged,metric.ls)

#M-H distance
diffDist5Zscore.avged <- distZ_avger(diffDist5Zscore)

overallTest(diffDist5Zscore.avged,"distanceZ",zscore=T,tails=2)
t.tester(diffDist5Zscore.avged,'distanceZ')

