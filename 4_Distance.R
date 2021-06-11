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
pb_download('nets_mix_clean10kS.RData',
            dest="data",
            tag="data.v.1")
pb_download('nets_mix_clean10kS.RData',
            dest="data",
            tag="data.v.1")
load('data/nets_mix_clean10kS.RData')
load('data/nets_mix_clean10kO.RData')



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
