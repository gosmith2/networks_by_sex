## 5: Sex differences vs species diet breadth

## calculates the rarefied degree of each species within the final 
## networks. Loads zscores and counts the number of significant 
## differences in network role for each observation.  Finally, the 
## number of differences is regressed against the rarefied degree.

library(lme4)
library(fossil)
library(tidyverse)

load('data/zscore50_5.RData')
load('data/spec_h.RData')
load('data/spec_y.RData')
load('data/spec_s.RData')



#calculate rarefied degree
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

save(rare.df,file='data/rare.RData')


#calculate number of differences, merge dataframes
metric.ls <- c("degree","weighted.betweenness",
               "weighted.closeness","d")



zjoin0 <- zscore_rare_joiner(zscore50_5.df,diffDist5Zscore,metric.ls)
zjoinST <- zscore_rare_joiner(zscore50_5bootST.cull,diffDist5ZscorebootST.cull,nets_mix_clean1kbootST,metric.ls)
zjoinST$Family <- zjoin0$Family[match(zjoinST$Genus,zjoin0$Genus)]
zjoinST$dataset <- 
zjoin5 <- zscore_rare_joiner(zscore50_5boot5.cull,diffDist5Zscoreboot5.cull,metric.ls)
zjoin5$Family <- zjoin0$Family[match(zjoin5$Genus,zjoin0$Genus)]
zjoin10<- zscore_rare_joiner(zscore50_5boot10.cull,diffDist5Zscoreboot10.cull,metric.ls)
zjoin10$Family <- zjoin0$Family[match(zjoin10$Genus,zjoin0$Genus)]
zjoin20<- zscore_rare_joiner(zscore50_5boot20.cull,diffDist5Zscoreboot20.cull,metric.ls)
zjoin20$Family <- zjoin0$Family[match(zjoin20$Genus,zjoin0$Genus)]


zjoin0$dist <- diffDist5Zscore$distanceZ[match(zjoin0$SpSiteYr,diffDist5Zscore$SpSiteYr)]

#model

#n0
rare.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zjoin0)
summary(rare.lmer)

rare.deg <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare.deg)

rare.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare.d)

rare.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare.clo)

rare.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare.md)

#more differences (across all metrics except MD) in more generalized species
#more generalized species had lower degree zscores
#more generalized species didn't differ in d'
#more generalized species had higher closeness zscores
#more generalized species had higher md zscores
#all pretty weak seeming

#ST
rareST.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zjoinST)
summary(rareST.lmer)

rareST.deg <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoinST,na.action=na.omit)
summary(rareST.deg)

rareST.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoinST,na.action=na.omit)
summary(rareST.d)

rareST.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr+~1|Family,data=zjoinST,na.action=na.omit)
summary(rareST.clo)

rareST.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoinST,na.action=na.omit)
summary(rareST.md)

zjoinST.melt <- melt(zjoinST,id.vars=c('SpSiteYr',"GenusSpecies","Genus","SiteYr","rareDeg","trials","plants"))

lineplot.CI(response=value,x.factor=rareDeg,data=zjoinST.melt,group=variable)


#n5
rare5.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zjoin5)
summary(rare5.lmer)

rare5.deg <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.deg)

rare5.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.d)

rare5.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.clo)

rare5.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.md)

#num of diffs had singularity issue
#more generalized species had lower degree zscores
#more generalized species had higher d'
#more generalized species had higher closeness zscores
#more generalized species had higher md zscores
#stronger than n0


#n10
rare10.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zjoin10)
summary(rare10.lmer)

rare10.deg <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.deg)

rare10.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.d)

rare10.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.clo)

rare10.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.md)

#num of diffs had singularity issue
#more generalized species had lower degree zscores
#more generalized species had higher d'
#more generalized species had higher closeness zscores
#more generalized species had higher md zscores
#stronger than n5


#n20
rare20.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zjoin20)
summary(rare20.lmer)

rare20.deg <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.deg)

rare20.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.d)

rare20.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.clo)

rare20.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.md)

#num of diffs had singularity issue
#more generalized species had lower degree zscores
#more generalized species had higher d'
#more generalized species had higher closeness zscores
#more generalized species had higher md zscores
#stronger
#adding plants as an interaction term doesnt do much




#----------------------
#family
#n0
rare0.degF <- lme(degree~Family,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare0.degF)
#halictids slightly less
rare0.dF <- lme(d~Family,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare0.dF)
#nothing
rare0.cloF <- lme(weighted.closeness~Family,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare0.cloF)
#halictids slightly more
rare0.mdF <- lme(dist~Family,random = ~1|SiteYr,data=zjoin0,na.action=na.omit)
summary(rare0.mdF)
#nothing

#n5
rare5.degF <- lme(degree~Family,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.degF)
#halictids less
rare5.dF <- lme(d~Family,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.dF)
#nothing
rare5.cloF <- lme(weighted.closeness~Family,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.cloF)
#halictids more
rare5.mdF <- lme(dist~Family,random = ~1|SiteYr,data=zjoin5,na.action=na.omit)
summary(rare5.mdF)
#halictids more, megachilids less


#n10
rare10.degF <- lme(degree~Family,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.degF)
#halictids less
rare10.dF <- lme(d~Family,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.dF)
#megachilids marg less
rare10.cloF <- lme(weighted.closeness~Family,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.cloF)
#halictids more
rare10.mdF <- lme(dist~Family,random = ~1|SiteYr,data=zjoin10,na.action=na.omit)
summary(rare10.mdF)
#halictids more



#n20
rare20.degF <- lme(degree~Family,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.degF)
#halictids less, megachilids less
rare20.dF <- lme(d~Family,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.dF)
#megachilids less
rare20.cloF <- lme(weighted.closeness~Family,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.cloF)
#halictids more
rare20.mdF <- lme(dist~Family,random = ~1|SiteYr,data=zjoin20,na.action=na.omit)
summary(rare20.mdF)
#halictids more, megachilids less
