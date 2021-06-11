## 5: Sex differences vs species diet breadth

## calculates the rarefied degree of each species within the final 
## networks. Loads and combines the node level zscores and tests 
## whether their distirbutions differ from 0. Finally, the 
## z-score values are regressed against the rarefied degree.

library(lme4)
library(fossil)
library(tidyverse)

load('data/zscore50_O.RData')
load('data/zscore50_S.RData')
load('data/diffDist5ZscoreO.RData')
load('data/diffDist5ZscoreS.RData')
load('data/spec_all.RData')
load('data/spec_h.RData')
load('data/spec_y.RData')
load('data/spec_s.RData')
load('data/nets_mix_clean10kO.RData')
load('data/nets_mix_clean10kS.RData')



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
load('data/rare.RData')

#Merge dataframes
metric2.ls <- c("degree","weighted.closeness","d")


zjoinO <- zscore_rare_joiner(zscore50_O.df,diffDist5ZscoreO,nets_mix_clean10kO,spec_all,metric2.ls,"O")
zjoinS <- zscore_rare_joiner(zscore50_S.df,diffDist5ZscoreS,nets_mix_clean10kS,spec_all,metric2.ls,"S")
#zjoinNEWO$dataset <- spec_all$dataset[match(zjoinNEWO$SiteYr,gsub(" ","_",spec_all$SiteYr))]
#zjoinNEWO$Family <- spec_all$Family[match(zjoinNEWO$GenusSpecies,spec_all$GenusSpecies)]
#zjoinNEWO$model <- rep("O")

###------------------
##Test: do the z-score distributions of the node-level metrics differ from zero?
  #no significant effect of dataset on any metric. all 4 metrics differed from 0 (intercept)


deg.O <- lme(degree~dataset,data=zjoinO,random=~1|SiteYr)
summary(deg.O)
clo.O <- lme(weighted.closeness~dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(clo.O)
d.O <- lme(d~dataset,data=zjoinO,random=~1|SiteYr)
summary(d.O)
md.O <- lme(dist~dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(md.O)


##genera where a difference in a network metric was observed
#unique(zjoinO$Genus[abs(zjoinO$degree)>1.96]) #9 genera
#unique(zjoinO$Genus[abs(zjoinO$d)>1.96]) #7 genera
#unique(zjoinO$Genus[abs(zjoinO$weighted.closeness)>1.96]) #10 genera
#unique(zjoinO$Genus[abs(zjoinO$dist)>1.96]) #10 genera


## --- validation with "S" dataset
  #no significant differences from 0 (intercept) or between datasets.
deg.S <- lme(degree~dataset,data=zjoinS,random=~1|SiteYr)
summary(deg.S)
clo.S <- lme(weighted.closeness~dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(clo.S)
d.S <- lme(d~dataset,data=zjoinS,random=~1|SiteYr)
summary(d.S)
md.S <- lme(dist~dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(md.S)




###------------------
##Test: Do males and females differ more in species that are more generalized 
## across the dataset (i.e., higher rarefied degree)?

deg.allp <- lme(degree ~ rareDeg+rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(deg.allp)
  #perfect. rare deg had neg? effect on degree (0.0256)
  #same neg effect with no plant#, everything else nothing

d.allp <- lme(d ~ rareDeg+rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(d.allp)

summary(d.allp)
  #no effects here
  #same, also no effects with plant

dist.allp <- lme(dist ~ rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)


summary(dist.allp)
#nothing here either
#without plant#, rareDeg sig (positive val, p 0.0002
#+plant -> sig positive

clo.all <- lme(weighted.closeness ~ rareDeg+plants+dataset,random=~1|SiteYr,data=zjoinO,na.action=na.omit)
summary(clo.all)
  #sig effect of dataset here
  #without plant#, closeness positive (p 0.0084), no effect of dataset
  #with plant+, sig of plants and sig of datasets

#-----
#S - no significant effects anywhere, regardless of plants or no plants
deg.allS <- lme(degree ~ rareDeg+rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(deg.allS)

d.allS <- lme(d ~ rareDeg+rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(d.allS)


dist.allS <- lme(dist ~ rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(dist.allS)

clo.allS <- lme(weighted.closeness ~ rareDeg+plants+dataset,random=~1|SiteYr,data=zjoinS,na.action=na.omit)
summary(clo.allS)




#---

deg.fam <- lme(dist ~ Family+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(deg.fam)


#zjoinST <- zscore_rare_joiner(zscore50_5bootST.cull,diffDist5ZscorebootST.cull,nets_mix_clean1kbootST,metric.ls)
#zjoinST$Family <- zjoin0$Family[match(zjoinST$Genus,zjoin0$Genus)]
#zjoinST$dataset <- 
#zjoin5 <- zscore_rare_joiner(zscore50_5boot5.cull,diffDist5Zscoreboot5.cull,metric.ls)
#zjoin5$Family <- zjoin0$Family[match(zjoin5$Genus,zjoin0$Genus)]
#zjoin10<- zscore_rare_joiner(zscore50_5boot10.cull,diffDist5Zscoreboot10.cull,metric.ls)
#zjoin10$Family <- zjoin0$Family[match(zjoin10$Genus,zjoin0$Genus)]
#zjoin20<- zscore_rare_joiner(zscore50_5boot20.cull,diffDist5Zscoreboot20.cull,metric.ls)
#zjoin20$Family <- zjoin0$Family[match(zjoin20$Genus,zjoin0$Genus)]


zjoin0$dist <- diffDist5Zscore$distanceZ[match(zjoin0$SpSiteYr,diffDist5Zscore$SpSiteYr)]

#model


#NEW

rare.deg <- lme(degree~rareDeg+model,random = ~1|SiteYr,data=zjoinNEW_all,na.action=na.omit)
summary(rare.deg)
#more generalized species had lower degree zscores

rare.d <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoinNEWO,na.action=na.omit)
summary(rare.d)
#more generalized species didn't differ in d'

rare.clo <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoinNEWO,na.action=na.omit)
summary(rare.clo)
#more generalized species had higher closeness zscores

rare.md <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoinNEWO,na.action=na.omit)
summary(rare.md)
#more generalized species had higher md zscores
#all pretty weak seeming


attach(zjoinNEWO);plot(plants,degree,col=c('red','blue')[as.factor(dataset)]);detach(zjoinNEWO)



#---S
rare.degS <- lme(degree~rareDeg,random = ~1|SiteYr,data=zjoinNEWS,na.action=na.omit)
summary(rare.degS)
#nuffin
rare.dS <- lme(d~rareDeg,random = ~1|SiteYr,data=zjoinNEWS,na.action=na.omit)
summary(rare.dS)
#nuffin
rare.cloS <- lme(weighted.closeness~rareDeg,random = ~1|SiteYr,data=zjoinNEWS,na.action=na.omit)
summary(rare.cloS)
#nuffin
rare.mdS <- lme(dist~rareDeg,random = ~1|SiteYr,data=zjoinNEWS,na.action=na.omit)
summary(rare.mdS)
#nuffin



rare.degO <- lme(degree~rareDeg+dataset+Family,random = ~1|SiteYr,data=zjoinNEWO,na.action=na.omit)
summary(rare.degO)

rare20.degF <- lme(degree~Family,random = ~1|SiteYr,data=zjoinNEWO,na.action=na.omit)
summary(rare20.degF)





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
