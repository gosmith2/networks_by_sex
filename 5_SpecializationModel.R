## 5: Sex differences vs species diet breadth

## calculates the rarefied degree of each species within the final 
## networks. Loads and combines the node level zscores and tests 
## whether their distirbutions differ from 0. Finally, the 
## z-score values are regressed against the rarefied degree.

library(nlme)
library(fossil)
library(tidyverse)
library(stringr)
source('prepNets.R')


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

save(zjoinO,file='data/zjoinO.RData')
save(zjoinS,file='data/zjoinS.RData')


#pb_upload('data/zjoinO.RData',
#          name="zjoinO.RData",
#          tag='data.v.1')
#pb_upload('data/zjoinS.RData',
#          name="zjoinS.RData",
#          tag='data.v.1')

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

deg.allO <- lme(degree ~ rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(deg.allO)

d.allO <- lme(d ~ rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(d.allO)

dist.allO <- lme(dist ~ rareDeg+plants+dataset,data=zjoinO,random=~1|SiteYr,na.action=na.omit)
summary(dist.allO)

clo.allO <- lme(weighted.closeness ~ rareDeg+plants+dataset,random=~1|SiteYr,data=zjoinO,na.action=na.omit)
summary(clo.allO)


#-----
#S - no significant effects anywhere
deg.allS <- lme(degree ~ rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(deg.allS)

d.allS <- lme(d ~ rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(d.allS)

dist.allS <- lme(dist ~ rareDeg+plants+dataset,data=zjoinS,random=~1|SiteYr,na.action=na.omit)
summary(dist.allS)

clo.allS <- lme(weighted.closeness ~ rareDeg+plants+dataset,random=~1|SiteYr,data=zjoinS,na.action=na.omit)
summary(clo.allS)

