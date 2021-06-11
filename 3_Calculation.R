#3: Calculations and testing

#Runs analyses on the species-level network parameters calculated in
#2_netBuilding.R to compare males and females within the observed and
#simulated networks. 

################################################

#Species level network traits

################################################

library(parallel)
library(piggyback)
library(bipartite)
library(tidyverse)
source('prepNets.R')

####-------------------------------------------####
#### Calculate node-level network parameters
####-------------------------------------------####

## calculate network stats at the node level, output into usable data frame
# reshuffling threshold is 5 in this case (i.e. there must be at least 5
# males and 5 females of a given species in a given network for it to be included)

metric.ls <- c("degree","closeness","d")

sex_trts_mix5O <- mclapply(nets_mix_clean10kO,
                           function(x) calcSpec(x, indiv = 5,index=metric.ls),
                           mc.cores = cores)

sex_trts_mix5S <- mclapply(nets_mix_clean10kS,
                           function(x) calcSpec(x, indiv = 5,index=metric.ls),
                           mc.cores = cores)



## confirm that the values are different
ifelse(any(sex_trts_mix5O[[1]]$weighted.closeness!=
             sex_trts_mix5O[[2]]$weighted.closeness),
       print("SUCCESS: the network statistics are different between randomizations"),
       print("WARNING: the network statistics are not different between randomizations")
)

#save and upload
save(sex_trts_mix5O,file='data/sex_trts_mix5O.RData')
save(sex_trts_mix5S,file='data/sex_trts_mix5S.RData')

pb_upload("data/sex_trts_mix5O.RData",
          name="sex_trts_mix5O.RData",
          tag="data.v.1")

##specify the metrics I'll be looking at, number of cores to use
cores <- 10

##clean up the dataset, add some necessary columns to speed things up
traits5S.ls <-
  mclapply(sex_trts_mix5S, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  return(x)
},mc.cores=cores)

traits5O.list <-
  mclapply(sex_trts_mix5O, function(x){
    x$SiteYr <- paste(x$Site, x$Year, sep="_")
    x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
    x <- filter(x,x$sex == "m" | x$sex == "f")
    return(x)
  },mc.cores=cores)


traits5S.list <- mclapply(c(1:10001),function(x){
  traits5S.ls[[x]]$weighted.closeness <- traits5S.lsc[[x]]$weighted.closeness
  return(traits5S.ls[[x]])
},mc.cores=cores)



###------------------
## calculate how different males and females are within each species
## within each SiteYr in each iteration. Comparison is 
## [male value - female value]: in the sexDiffs output, large values
## indicate that males had larger values than females.

metric2.ls <- c('degree','weighted.closeness','d')

sexDiffs5O.df <- makeComp(traits5O.list, metric2.ls, comparison = "diff")
sexDiffs5S.df <- makeComp(traits5S.list, metric2.ls, comparison = "diff")


## Saving and uploading after that long step
save(sexDiffs5O.df, file = 'data/sexDiffs5O.RData')
save(sexDiffs5S.df, file = 'data/sexDiffs5S.RData')


pb_upload('data/sexDiffs5S.RData',
          name='sexDiffs5S.RData',
          tag="data.v.1")

pb_upload('data/sexDiffs5O.RData',
          name='sexDiffs5O.RData',
          tag="data.v.1")


## **************************************************************** 
##Calculate how different the observed values were from the simulated values.
## when zscore=F, outupts the proportion of simulations where the value
## was less than or 50% equal to the obesrved difference value. 
## High proportions indicate that males were larger than females to a
## greater degree than the null expectation. Low proportions indicate that 
## females had higher values than males to a greater degree than expected

## downloading the output to skip above step 
#pb_download("sexDiffs5O.df",
#            dest="data",
#            tag="data.v.1")
#load('data/sexDiffs5boot.df')

#sexDiffsProp50_5.df <- calcNullProp50(sexDiffs5.df,
#                                      metric.ls,
#                                      zscore=FALSE)

##generate z-scores
zscore50_O.df <- calcNullProp50(sexDiffs5O.df, 
                                metric2.ls,
                                zscore=TRUE)
zscore50_S.df <- calcNullProp50(sexDiffs5S.df, 
                                   metric2.ls,
                                   zscore=TRUE)

save(zscore50_O.df,file='data/zscore50_O.RData')
save(zscore50_S.df,file='data/zscore50_S.RData')


pb_upload('data/zscore50_S.RData',name='zscore50_S.RData',tag='data.v.1')
pb_upload('data/zscore50_O.RData',name='zscore50_O.RData',tag='data.v.1')



pb_download('zscore50_O.RData',dest='data',tag='data.v.1')








#genera where a difference in a network metric was observed
unique(zscore50_5O.df$Genus[abs(zscore50_5.df$degree)>1.96])
unique(zscore50_5O.df$Genus[abs(zscore50_5.df$d)>1.96])
unique(zscore50_5O.df$Genus[abs(zscore50_5.df$weighted.closeness)>1.96])


######-------------------------------

#save and upload
save(zscore50_5.df,file="data/zscore50_5.RData")
#save(sexDiffsProp50_5.df,file='data/sexDiffsProp50_5.Rdata')

pb_upload("data/zscore50_5bootST.RData",
          name="zscore50_5bootST.RData",
          tag="data.v.1")




t1 <- lme(d~dataset,data=zjoinNEWO,random=~1|SiteYr,na.action=na.omit)
summary(t1)

t2 <- lme.maker(zjoinNEWO,metric.ls,random.eff='SiteYr')

#bootstuff

#
zscore50_5NEWOST.cull <- subset(zscore50_NEWOST.df,zscore50_NEWOST.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5NEWSST.cull <- subset(zscore50_NEWSST.df,zscore50_NEWSST.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5NEWDST.cull <- subset(zscore50_NEWDST.df,zscore50_NEWDST.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)

diffDist5ZscoreNEWSST.cull <- subset(diffDist5ZscoreNEWSST,gsub("\\.","_",diffDist5ZscoreNEWSST$Level) %in% diffDist5Zscore$Level)
diffDist5ZscoreNEWDST.cull <- subset(diffDist5ZscoreNEWDST,gsub("\\.","_",diffDist5ZscoreNEWDST$Level) %in% diffDist5Zscore$Level)
diffDist5ZscoreNEWOST.cull <- subset(diffDist5ZscoreNEWOST,gsub("\\.","_",diffDist5ZscoreNEWOST$Level) %in% diffDist5Zscore$Level)


#
zscore50_5bootST.cull <- subset(zscore50_5bootST.df,zscore50_5bootST.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5bootT.cull <- subset(zscore50_5bootT.df,zscore50_5bootT.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5boot5.cull <- subset(zscore50_5boot5.df,zscore50_5boot5.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5boot10.cull <- subset(zscore50_5boot10.df,zscore50_5boot5.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)
zscore50_5boot20.cull <- subset(zscore50_5boot20.df,zscore50_5boot5.df$SpSiteYr %in% zscore50_5.df$SpSiteYr)

nodeboot <- boot_sumarizer(c(1:5),list(zscore50_5.df,zscore50_5bootST.cull,zscore50_5boot5.cull,zscore50_5boot10.cull,zscore50_5boot20.cull),c("Zero",'ST','Five','Ten','Twenty'),metric.ls)

diffDist5ZscorebootSTJ.cull <- subset(diffDist5ZscorebootSTJ,gsub("\\.","_",diffDist5ZscorebootSTJ$Level) %in% diffDist5ZscoreJ$Level)
diffDist5ZscorebootT.cull <- subset(diffDist5ZscorebootT,gsub("\\.","_",diffDist5ZscorebootT$Level) %in% diffDist5Zscore$Level)
diffDist5Zscoreboot5J.cull <- subset(diffDist5Zscoreboot5J,gsub("\\.","_",diffDist5Zscoreboot5J$Level) %in% diffDist5ZscoreJ$Level)
diffDist5Zscoreboot10J.cull <- subset(diffDist5Zscoreboot10J,gsub("\\.","_",diffDist5Zscoreboot10J$Level) %in% diffDist5ZscoreJ$Level)
diffDist5Zscoreboot20J.cull <- subset(diffDist5Zscoreboot20J,gsub("\\.","_",diffDist5Zscoreboot20J$Level) %in% diffDist5ZscoreJ$Level)

distboot <- boot_sumarizer(c(1:5),list(diffDist5ZscoreJ,diffDist5ZscorebootSTJ.cull,diffDist5Zscoreboot5J.cull,diffDist5Zscoreboot10J.cull,diffDist5Zscoreboot20J.cull),c("Zero",'ST','Five','Ten','Twenty'),'distanceZ')

sexlvlProp50ZbootTI.cull <- subset(sexlvlProp50ZbootTI,sexlvlProp50ZbootTI$SpSiteYr %in% sexlvlProp50Z$SpSiteYr)
sexlvlProp50ZbootST.cull <- subset(sexlvlProp50ZbootST,sexlvlProp50ZbootST$SpSiteYr %in% sexlvlProp50Z$SpSiteYr)
sexlvlProp50Zboot5.cull <- subset(sexlvlProp50Zboot5,sexlvlProp50Zboot5$SpSiteYr %in% sexlvlProp50Z$SpSiteYr)
sexlvlProp50Zboot10.cull <- subset(sexlvlProp50Zboot,sexlvlProp50Zboot$SpSiteYr %in% sexlvlProp50Z$SpSiteYr)
sexlvlProp50Zboot20.cull <- subset(sexlvlProp50Zboot20,sexlvlProp50Zboot20$SpSiteYr %in% sexlvlProp50Z$SpSiteYr)

sexmet <- colnames(sexlvlProp50Z)[1:8]
netboot <- boot_sumarizer(c(1:5),list(sexlvlProp50Z,sexlvlProp50ZbootST.cull,sexlvlProp50Zboot5.cull,sexlvlProp50Zboot10.cull,sexlvlProp50Zboot20.cull),c("Zero",'ST','Five','Ten','Twenty'),sexmet)

all.boot <- rbind(nodeboot,distboot,netboot)
all.boot$label <- factor(all.boot$label,levels = c('Zero','ST',"Five",'Ten','Twenty'))

lineplot.CI(x.factor=label,response=mean,group=metric,data=all.boot)
