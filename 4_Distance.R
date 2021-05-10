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

load("data/nets_mix_clean1kbootST.Rdata")

#calculate the M-H index for species with at least 5 males and 5 females
distValues5J.df <- distComp(nets_mix_clean2k[1:2],"jaccard",indiv=5) #try chao, try jaccard<-

#add siteyr column, remove NAs from sites that did not have species with the 
#minimum abundances
distValues5J.df$SpSiteYr <- paste(distValues5J.df$GenusSpecies,
                                 gsub("\\.","_",distValues5J.df$SiteYr), sep="_")
distValues5J.df <- filter(distValues5J.df,distance!='NA')


#Save and upload distances
save(distValues5.df,file='data/distValues5.RData')

pb_upload("data/distValues5boot20.RData",
          name="distValues5boot20.RData",
          tag="data.v.1")

## ****************************************************************
#download to the previous skip step
pb_download("distValues5.RData",
            dest="data",
            tag="data.v.1")

#compare observed M-H distances to simulated null networks, generating zscores
diffDist5ZscoreJ <- calcDistZ(distValues5J.df,"SpSiteYr",zscore=T)
save(diffDist5ZscorebootT,file='data/diffDist5ZscorebootT.RData')
pb_upload("data/diffDist5ZscorebootT.RData",
          name="diffDist5ZscorebootT.RData",
          tag="data.v.1")
#diffDist5 <- calcDistZ(distValues5.df,"SpSiteYr",zscore=F)


## Test: proportion of species+sites where M-H distance in
## observed network was different than many of the simulations. The output
## is the proportion of observations (Sp+Site+Year) where the distance
## diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

#overallTest(diffDist5,"distanceZ",zscore=F,tails=2)
overallTest(diffDist5Zscoreboot20,"distanceZ",zscore=T,tails=2)

t.tester(diffDist5Zscoreboot20,'distanceZ')
  #mean for distribution of zscore values diff from 0 
  
#genera where a difference in M-H was observed
unique(word(unique(diffDist5ZscorebootST$Level[abs(diffDist5Zscore$distanceZ)>1.96]),1))
  #10 genera
  #with ST, there are 27 genera
  
pb_upload("data/diffDist5ZscorebootST.RData",
          name="diffDist5ZscorebootST.RData",
          tag="data.v.1")



######-------Section S1: average together sites within a given dataset

load('data/zscore50_5.RData')
load('data/diffDist5Zscore.RData')
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

