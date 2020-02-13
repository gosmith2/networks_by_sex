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

overallTest(sexDiffsProp50_2.df, metric.ls, zscore=F)

overallTest(zscore50_2.df, metric.ls, zscore=T)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=2)

#splevel
spLevelTest(zscore50_2.df,metric.ls,zscore=T)
#results: tons of zeros

#familylevel
spLevelTest(zscore50_2.df,metric.ls,zscore=T,level="Family")


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

plot(density(zscore50_2.df$weighted.closeness, na.rm = T))
abline(v=0)
#regress amt of difference against absolute specialization?
#i.e, are more generalized sp more different b/w males and females?


#################################################

#model for divergence vs generality

#################################################

load("data/zscore50_2.RData")
polltraitsSF.df<-read.csv("data/polltraitsSF.csv")
polltraitsY.df<-read.csv("data/polltraitsY.csv")

zscore50_2SI.df %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  zscore50_2SI.df

zscore_traits.df <- inner_join(zscore50_2.df, polltraitsSF.df, by="GenusSpecies")

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



##
#----------Roswell paper notes

#generally, a lot more sampling, I think. like, all day
  #for 3 days

#Do males and females overlap in diet?
  #9999 iterations: this many reqed to stabilize for 0.05 alpha
    #based on North, Curtis, and Sham 2002
  #)"When the observed dissimilarity was greater than 9500
    #ofthe 9999 simulated dissimilarities, we concluded that
    #we had detected a difference in the pattern of floral
    #visitation between conspecific male and female bees, 
    #given the observed diet breadth and abundance of each sex.
  #dissimilarity calculated using: Morisita-Horn index
    #good for mixed-size and small sample data
  #comparison at the SPECIES level, across all sites and 
    #sample rounds (several w/in 2016)




