#6: Modeling based on traits

library(fossil)
load("data/zscore50_2.RData")
polltraitsSF.df<-read.csv("data/polltraitsSF.csv")
polltraitsY.df<-read.csv("data/polltraitsY.csv")

#make COMPLETE network for each dataset as a trainer

wholeNet <- function(sampData){
  obs <- sampData[sampData$Sex=="f",]
  agg.spec <- aggregate(list(abund=obs$GenusSpecies),
                         list(GenusSpecies=
                                obs$GenusSpecies,
                              PlantGenusSpecies=
                                obs$PlantGenusSpecies),
                         length)
  net <- samp2site.spp(agg.spec$PlantGenusSpecies, 
                       agg.spec$GenusSpecies,
                       agg.spec$abund)
  return(net)
}

calcRareDeg <- function(sampData, name){
  net <- wholeNet(sampData)
  rare <- apply(net,2,chao1)
  rare.df <- data.frame(GenusSpecies = names(rare),
                        r.deg = rare,
                        row.names = NULL)
  names(rare.df) <- c("GenusSpecies",name)
  return(rare.df)
}

rare.hr <- calcRareDeg(spec.h,"rare.hr")
rare.si <- calcRareDeg(spec.s,"rare.si")
rare.yo <- calcRareDeg(spec.y,"rare.yo")

rare.df <- full_join(rare.hr,rare.si,by="GenusSpecies")
rare.df <- full_join(rare.df,rare.yo,by="GenusSpecies")
rare.df$r.deg1 <- ifelse(is.na(rare.df$rare.hr),
                        ifelse(is.na(rare.df$rare.si),
                               rare.df$rare.yo,
                               rare.df$rare.si),
                        rare.df$rare.hr)





#zscore50_2.df %>%
#  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
#  zscore50_2.df

#zscore_traits.df <- inner_join(zscore50_2.df, polltraitsSF.df, by="GenusSpecies")
#
#zscore_traits.df$distance <- distValues.df$distance[match(zscore_traits.df$SpSiteYr,
#                                                          distValues.df$SpSiteYr)]


traitFrame <- function(diffs, dist) {
  #extract species from SpSiteYr
  diffs %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  diff1
  
  #add in the rarefied degrees, matching by species
  diff1 <- left_join(diff1, rare.df, by="GenusSpecies")
  
  #add in distance index proportions, matching by SpSiteYr
  diff1$distance <- dist$distanceZ[match(diff1$SpSiteYr,
                                      dist$Level)]
  res <- str_match(diff1$SpSiteYr, "_(.*?)_")
  diff1$Site <- res[,2]
  diff1$r.deg <- ifelse(diff1$Site %in% unique(spec.h$Site), 
                        diff1$rare.hr,
                        ifelse(diff1$Site %in% unique(spec.s$Site),
                               diff1$rare.si,
                               diff1$rare.yo)
  )
return(diff1[,-c(8:11)])
}

library(stringr)
res <- str_match(sexDiffsProp50_2.df$SpSiteYr, "_(.*?)_")
res[,2]

diff_traits2.df<-traitFrame(sexDiffsProp50_2.df,diffDist2)
diff_traits5.df<-traitFrame(sexDiffsProp50_5.df,diffDist5)
diff_traits10.df<-traitFrame(sexDiffsProp50_10.df,diffDist10)

diff_traits5.df$degree.obs <- sexDiffs5.df[[1]]$degree
diff_traits5.df$str.obs <- sexDiffs5.df[[1]]$species.strength
diff_traits5.df$btw.obs <- sexDiffs5.df[[1]]$weighted.betweenness
diff_traits5.df$clo.obs <- sexDiffs5.df[[1]]$weighted.closeness
diff_traits5.df$d.obs <- sexDiffs5.df[[1]]$d
diff_traits5.df$dist.obs <- distValues5obs.df$distance[match(diff_traits5.df$SpSiteYr,
                                                             distValues5obs.df$SpSiteYr)]



##########################################

#diff thresholds and specialization models

##########################################

###########
#r.deg vs degree: 2 and 5 sig, 10 marg
rareDeg2<-lme(degree~r.deg,data=diff_traits2.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(rareDeg2)
plot(diff_traits2.df$degree~diff_traits2.df$r.deg)
abline(lm(diff_traits2.df$degree~diff_traits2.df$r.deg))


rareDeg5<-lme(degree~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareDeg5)
plot(diff_traits5.df$degree~diff_traits5.df$r.deg)
abline(lm(diff_traits5.df$degree~diff_traits5.df$r.deg))

rareDeg10<-lme(degree~r.deg,data=diff_traits10.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareDeg10)
plot(diff_traits10.df$degree~diff_traits10.df$r.deg)


rareDeg5obs<-lme(degree.obs~r.deg,data=diff_traits5.df,
                            random=~1|GenusSpecies,na.action=na.omit)
summary(rareDeg5obs)
plot(diff_traits5.df$degree.obs~diff_traits5.df$r.deg)
abline(lm(diff_traits5.df$degree.obs~diff_traits5.df$r.deg))

###############
#r.deg vs str: all 3 sig
rareStr2<-lme(species.strength~r.deg,data=diff_traits2.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareStr2)

rareStr5<-lme(species.strength~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareStr5)
plot(diff_traits5.df$species.strength~diff_traits5.df$r.deg)
abline(lm(diff_traits5.df$species.strength~diff_traits5.df$r.deg))

rareStr10<-lme(species.strength~r.deg,data=diff_traits10.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareStr10)
plot(diff_traits10.df$degree~diff_traits10.df$r.deg)


rareStr5obs<-lme(str.obs~r.deg,data=diff_traits5.df,
                 random=~1|GenusSpecies,na.action=na.omit)
summary(rareStr5obs)
plot(diff_traits5.df$str.obs~diff_traits5.df$r.deg)
abline(lm(diff_traits5.df$str.obs~diff_traits5.df$r.deg))

###################
#r.deg vs b/w: none sig
rareBtw2<-lme(weighted.betweenness~r.deg,data=diff_traits2.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareBtw2)

rareBtw5<-lme(weighted.betweenness~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareBtw5)

rareBtw10<-lme(weighted.betweenness~r.deg,data=diff_traits10.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareBtw10)

#marg
rareBtw5obs<-lme(btw.obs~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareBtw5obs)

###################
#r.deg vs closeness: 2 sig, others not
rareClo2<-lme(weighted.closeness~r.deg,data=diff_traits2.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareClo2)

rareClo5<-lme(weighted.closeness~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareClo5)

rareClo10<-lme(weighted.closeness~r.deg,data=diff_traits10.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareClo10)


rareClo5obs<-lme(clo.obs~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareClo5obs)


###################
#r.deg vs d: none sig
rareD2<-lme(d~r.deg,data=diff_traits2.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareD2)

rareD5<-lme(d~r.deg,data=diff_traits5.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareD5)

rareD10<-lme(d~r.deg,data=diff_traits10.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareD10)

rareD5obs<-lme(d.obs~r.deg,data=diff_traits5.df,
            random=~1|GenusSpecies,na.action=na.omit)
summary(rareD5obs)


############################
#distance by r.deg: all three sig
rareDist2<-lme(distance~r.deg,data=diff_traits2.df,
                random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist2)
plot(distance~r.deg,data=diff_traits2.df)

rareDist5<-lme(distance~r.deg,data=diff_traits5.df,
                random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist5)
plot(distance~r.deg,data=diff_traits5.df)
abline(lm(distance~r.deg,data=diff_traits5.df))

rareDist10<-lme(distance~r.deg,data=diff_traits10.df,
                 random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist10)


rareDist5obs<-lme(dist.obs~r.deg,data=diff_traits5.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist5obs)



















#################################################

#model for divergence vs generality

#################################################

#Lecty models
lectyDeg<-lme(degree.x~Lecty,data=diff_traits2.df,
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

lectyDist<-lme(distance~Lecty,data=diff_traits2.df,
               random=~1|GenusSpecies,na.action=na.omit)
summary(lectyDist)
plot(diff_traits2.df$distance~diff_traits2.df$Lecty)

rareDist<-lme(distance~r.degree,data=diff_traits2.df,
              random=~1|GenusSpecies,na.action=na.omit)
summary(rareDist)
plot(diff_traits2.df$distance~diff_traits2.df$r.degree)

degDist<-lme(distance~degree.y,data=zscore_traits.df,
             random=~1|GenusSpecies,na.action=na.omit)
summary(degDist)
plot(zscore_traits.df$distance~zscore_traits.df$degree.y)

dDist<-lme(distance~d.y,data=zscore_traits.df,
           random=~1|GenusSpecies,na.action=na.omit)
summary(dDist)
plot(zscore_traits.df$distance~zscore_traits.df$d.y)



