## 5: Dissimilarity indexing

library(vegan)

load("data/mix_netsYHS.RData")

distValues.df <- distComp(nets.mix.clean,"horn")
distValues2.df <- distComp(nets.mix.clean,"horn",indiv=2)
distValues5.df <- distComp(nets.mix.clean,"horn",indiv=5)
distValues10.df <- distComp(nets.mix.clean,"horn",indiv=10)

distValues5obs.df<-subset(distValues5.df,distValues5.df$sim==1)

distValues10.df$SpSiteYr <- paste(distValues10.df$GenusSpecies,
                                distValues10.df$SiteYr,
                                sep="_")
distValues10.df$SpSiteYr <- gsub("\\.","_",distValues10.df$SpSiteYr)

save(distValues.df,file='data/distValues.RData')
save(distValues2.df,file='data/distValues2.RData')
save(distValues5.df,file='data/distValues5.RData')
save(distValues10.df,file='data/distValues10.RData')

pb_upload("data/distValues10.RData",
          name="distValues10.RData",
          tag="data.v.1")
pb_download("distValues10.RData",
            dest="data",
            tag="data.v.1")

zscoreDist <- calcDistZ(distValues.df,"SpSiteYr",zscore=T)
diffDist <- calcDistZ(distValues.df,"SpSiteYr",zscore=F)
diffDist2 <- calcDistZ(distValues2.df,"SpSiteYr",zscore=F)
diffDist5 <- calcDistZ(distValues5.df,"SpSiteYr",zscore=F)
diffDist10 <- calcDistZ(distValues10.df,"SpSiteYr",zscore=F)



overallTest(diffDist,"distanceZ",zscore=F,tails=-1)
overallTest(diffDist2,"distanceZ",zscore=F,tails=-1)
overallTest(diffDist5,"distanceZ",zscore=F,tails=-1)
overallTest(diffDist10,"distanceZ",zscore=F,tails=-1)



#hive plot


