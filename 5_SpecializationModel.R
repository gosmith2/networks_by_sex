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


#calculate number of differences, merge dataframes
metric.ls <- c("degree","weighted.betweenness",
               "weighted.closeness","d")

zscore50_5.df$GenusSpecies <- word(zscore50_5.df$SpSiteYr,1,sep = "_")
zscore50_5.df$SiteYr <- word(zscore50_5.df$SpSiteYr,2,3,sep='_')
zscore50_5.df$diffs <- unlist(lapply(zscore50_5.df$SpSiteYr,function(x){
  sp <- zscore50_5.df[zscore50_5.df$SpSiteYr==x,]
  val <- sum(sp[metric.ls]>1.96|sp[metric.ls]<(-1.96), na.rm=T)
  return(val)
  })
)
zscore50_5.df$rareDeg <- rare.df$r.deg1[match(zscore50_5.df$GenusSpecies,rare.df$GenusSpecies)]
zscore50_5.df$trials <- 4


#model
rare.lmer <- glmer(cbind(diffs,trials) ~ rareDeg +(1|SiteYr), family = binomial, data=zscore50_5.df)
summary(rare.lmer)
