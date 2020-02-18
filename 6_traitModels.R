#6: Modeling based on traits

load("data/zscore50_2.RData")
polltraitsSF.df<-read.csv("data/polltraitsSF.csv")
polltraitsY.df<-read.csv("data/polltraitsY.csv")

zscore50_2.df %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  zscore50_2.df

zscore_traits.df <- inner_join(zscore50_2.df, polltraitsSF.df, by="GenusSpecies")

zscore_traits.df$distance <- distValues.df$distance[match(zscore_traits.df$SpSiteYr,
                                                          distValues.df$SpSiteYr)]
#################################################

#model for divergence vs generality

#################################################

#Lecty models
lectyDeg<-lme(degree.x~Lecty,data=zscore_traitsSI.df,
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
abline(lm(species.strength~r.degree,data=zscore_traits.df))


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

##################################

#Distance

##################################

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



