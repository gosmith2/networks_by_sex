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



zscore_traits.df$distance <- distValues.df$distance[match(zscore_traits.df$SpSiteYr,
                                              distValues.df$SpSiteYr)]

lectyDist<-lme(distance~Lecty,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(lectyDist)
plot(zscore_traits.df$distance~zscore_traits.df$Lecty)

rareDist<-lme(distance~r.degree,data=zscore_traits.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist)
plot(zscore_traits.df$distance~zscore_traits.df$r.degree)

degDist<-lme(distance~degree.y,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(degDist)
plot(zscore_traits.df$distance~zscore_traits.df$degree.y)

dDist<-lme(distance~d.y,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(dDist)
plot(zscore_traits.df$distance~zscore_traits.df$d.y)


