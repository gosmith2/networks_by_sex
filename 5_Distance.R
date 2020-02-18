## 5: Dissimilarity indexing

library(vegan)

load("data/mix_netsYHS.RData")

distValues.df <- distComp(nets.mix.clean,"horn")
distValues.df$SpSiteYr <- paste(distValues.df$GenusSpecies,
                                distValues.df$SiteYr,
                                sep="_")
distValues.df$SpSiteYr <- gsub("\\.","_",distValues.df$SpSiteYr)

save(distValues.df,file='data/distValues.RData')

pb_upload("data/distValues.RData",
          name="distValues.RData",
          tag="data.v.1")
pb_download("distValues.RData",
            dest="data",
            tag="data.v.1")

zscoreDist <- calcDistZ(distValues.df,"SpSiteYr")

overallTest(zscoreDist,"distanceZ")



#hive plot


