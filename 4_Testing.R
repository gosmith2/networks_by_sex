## 4: Testing and models


pb_download("zscore50_2.RData",
            dest="data",
            tag="data.v.1")

pb_download("sexDiffsProp50YH_2.Rdata",
            dest="data",
            tag="data.v.1")

load("data/zscore50_2.RData")
load("data/sexDiffsProp50YH_2.Rdata")



###------------------
##Test: proportion of species+sites where m v f difference in
##observed network was larger than many of the simulations

#same test as above, but with 50% >=, rather than >=

overallTest(sexDiffsProp50_2.df, metric.ls, zscore=F)

overallTest(zscore50_2.df, metric.ls, zscore=T)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=2)

spLevelTest(zscore50_2.df,metric.ls,zscore=T)
#results: tons of zeros, species seem to fluctuate around 0.05

sexDiffsProp50_2.df %>%
  mutate(Sp = gsub( "_.*$", "", SpSiteYr))%>%
  arrange(desc(Sp)) -> sexDiffsProp50_2.df

zscore50_2.df %>%
  mutate(Sp = gsub( "_.*$", "", SpSiteYr))





##density plot for degree, absolute differences
plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)
#regress amt of difference against absolute specialization?
#i.e, are more generalized sp more different b/w males and females?


#################################################

#model for divergence vs generality

#################################################

load("data/zscore50_2.RData")
polltraitsSF.df<-read.csv("data/polltraitsSF.csv")
polltraitsY.df<-read.csv("data/polltraitsY.csv")

zscore50_2.df %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  zscore50_2.df

#sp.ls <- unique(zscore50_2.df$Sp)

#polltraits<-bind_rows(polltraitsY.df[c("GenusSpecies","Lecty")],
#          polltraitsSF.df[c("GenusSpecies","Lecty")])

#zscore_traits.df <- inner_join(zscore50_2.df, polltraits, by="GenusSpecies")
zscore_traits.df <- inner_join(zscore50_2.df, polltraitsSF.df, by="GenusSpecies")

#Lecty models
lectyDeg<-lme(degree~Lecty,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(lectyDeg)

lectyStr<-lme(species.strength~Lecty,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(lectyStr)

lectyBtw<-lme(weighted.betweenness~Lecty,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(lectyBtw)

lectyClo<-lme(weighted.closeness~Lecty,data=zscore_traits.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(lectyClo)

lectyD<-lme(d.x~Lecty,data=zscore_traits.df,
            random=~1|GenusSpecies,na.action=na.omit)
summary(lectyD)


#r.deg models: full dataset
rareDeg<-lme(degree.x~r.degree,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(rareDeg)
plot(zscore_traits.df$degree.x~zscore_traits.df$r.degree)
abline(lm(degree.x~r.degree,data=zscore_traits.df))

rareStr<-lme(species.strength~r.degree,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(rareStr)
plot(zscore_traits.df$species.strength~zscore_traits.df$r.degree)


rareBtw<-lme(weighted.betweenness~r.degree,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(rareBtw)
plot(zscore_traits.df$weighted.betweenness~zscore_traits.df$r.degree)


rareClo<-lme(weighted.closeness~r.degree,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(rareClo)
plot(zscore_traits.df$weighted.closeness~zscore_traits.df$r.degree)
abline(lm(weighted.closeness~r.degree,data=zscore_traits.df))


#r.deg models, excluding the hypergeneralists
#without h. tripartitus and L. incompletum, effecs disappear
z_traits_spec <- subset(zscore_traits.df,zscore_traits.df$r.degree<75)

rDegSpec<-lme(degree.x~r.degree,data=z_traits_spec,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rDegSpec)

rStrSpec<-lme(species.strength~r.degree,data=z_traits_spec,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rStrSpec)

rBtwSpec<-lme(weighted.betweenness~r.degree,data=z_traits_spec,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rBtwSpec)

rCloSpec<-lme(weighted.closeness~r.degree,data=z_traits_spec,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rCloSpec)

rDSpec<-lme(d.x~r.degree,data=z_traits_spec,
            random=~1|GenusSpecies,na.action=na.omit)
summary(rDSpec)

#d': with full data, same patterns as r.degree.
dDeg<-lme(degree.x~d,data=zscore_traits.df,
          random=~1|GenusSpecies,na.action=na.omit)
summary(dDeg)

dStr<-lme(species.strength~d,data=zscore_traits.df,
          random=~1|GenusSpecies,na.action=na.omit)
summary(dStr)

dBtw<-lme(weighted.betweenness~d,data=zscore_traits.df,
          random=~1|GenusSpecies,na.action=na.omit)
summary(dBtw)

dClo<-lme(weighted.closeness~d,data=zscore_traits.df,
          random=~1|GenusSpecies,na.action=na.omit)
summary(dClo)
