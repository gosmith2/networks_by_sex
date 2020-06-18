## 4: Dissimilarity indexing

## Loads networks and calculates the Morisita-Horn dissimilarity
## index for males and females within each species within each network.
## Lastly, the number of observations showing large differences from 
## the null is calculated. 

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

#distValues5obs.df<-subset(distValues5.df,distValues5.df$sim==1)

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

######-------------------------------

# - Alternative code that will be removed before pub

######-------------------------------

distValues.df <- distComp(nets.mix.clean,"horn")
distValues2.df <- distComp(nets.mix.clean,"horn",indiv=2)
distValues10.df <- distComp(nets.mix.clean,"horn",indiv=10)


save(distValues.df,file='data/distValues.RData')
save(distValues2.df,file='data/distValues2.RData')
save(distValues10.df,file='data/distValues10.RData')


distValues10.df$SpSiteYr <- paste(distValues10.df$GenusSpecies,
                                  distValues10.df$SiteYr,
                                  sep="_")
distValues10.df$SpSiteYr <- gsub("\\.","_",distValues10.df$SpSiteYr)


zscoreDist <- calcDistZ(distValues.df,"SpSiteYr",zscore=T)
diffDist <- calcDistZ(distValues.df,"SpSiteYr",zscore=F)
diffDist2 <- calcDistZ(distValues2.df,"SpSiteYr",zscore=F)
diffDist10 <- calcDistZ(distValues10.df,"SpSiteYr",zscore=F)
overallTest(diffDist10,"distanceZ",zscore=F,tails=-1)
overallTest(diffDist,"distanceZ",zscore=F,tails=-1)
overallTest(diffDist2,"distanceZ",zscore=F,tails=-1)