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
pb_download('nets_mix_clean.RData',
            dest="data",
            tag="data.v.1")

load("data/nets_mix_clean.Rdata")

#calculate the M-H index for species with at least 5 males and 5 females
distValues5.df <- distComp(nets_mix_clean,"horn",indiv=5)

#add siteyr column, remove NAs from sites that did not have species with the 
#minimum abundances
distValues5.df$SpSiteYr <- paste(distValues5.df$GenusSpecies,
                                 distValues5.df$SiteYr)
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

#compare observed M-H distances to simulated null networks
diffDist5 <- calcDistZ(distValues5.df,"SpSiteYr",zscore=F)

#calculate again to generate z-scores for plotting
diffDist5Zscore <- calcDistZ(distValues5.df,"SpSiteYr",zscore=T)
save(diffDist5Zscore,file='data/diffDist5Zscore.RData')

## Test: proportion of species+sites where M-H distance in
## observed network was different than many of the simulations. The output
## is the proportion of observations (Sp+Site+Year) where the distance
## diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

overallTest(diffDist5,"distanceZ",zscore=F,tails=2)


pb_upload("data/diffDist5Zscore.RData",
          name="diffDist5Zscore.RData",
          tag="data.v.1")
