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

###------------------
## Setup
pb_download("sex_trts_mix5.RData",
            dest="data",
            tag="data.v.1")

load("data/sex_trts_mix5.RData")

##specify the metrics I'll be looking at, number of cores to use
metric.ls <- c("degree","weighted.betweenness",
               "weighted.closeness","d")
cores <- 10

##clean up the dataset, add some necessary columns to speed things up
traits5bootT.ls <-
  mclapply(sex_trts_mix5bootT, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  return(x)
},mc.cores=cores)


###------------------
## calculate how different males and females are within each species
## within each SiteYr in each iteration. Comparison is 
## [male value - female value]: in the sexDiffs output, large values
## indicate that males had larger values than females.

sexDiffs5bootT.df <- makeComp(traits5bootT.ls, metric.ls, comparison = "diff")

## Saving and uploading after that long step
save(sexDiffs5bootST.df, file = 'data/sexDiffs5bootST.RData')

pb_upload('data/sexDiffs5boot.RData',
          name='sexDiffs5boot.RData',
          tag="data.v.1")


## **************************************************************** 
##Calculate how different the observed values were from the simulated values.
## when zscore=F, outupts the proportion of simulations where the value
## was less than or 50% equal to the obesrved difference value. 
## High proportions indicate that males were larger than females to a
## greater degree than the null expectation. Low proportions indicate that 
## females had higher values than males to a greater degree than expected

## downloading the output to skip above step 
pb_download("sexDiffs5boot.df",
            dest="data",
            tag="data.v.1")
load('data/sexDiffs5boot.df')

#sexDiffsProp50_5.df <- calcNullProp50(sexDiffs5.df,
#                                      metric.ls,
#                                      zscore=FALSE)

##generate z-scores
zscore50_5bootT.df <- calcNullProp50(sexDiffs5bootT.df, 
                                metric.ls,
                                zscore=TRUE)
save(zscore50_5bootT.df,file='data/zscore50_5bootT.RData')
pb_upload('data/zscore50_5bootT.RData',name='zscore50_5bootT.RData',tag='data.v.1')

###------------------
##Test: proportion of species+sites where m v f difference in
##observed network was larger than many of the simulations. The output
##is the proportion of observations (Sp+Site+Year) where the male-female
## difference diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

#overallTest(sexDiffsProp50_5.df, metric.ls, tails=2, zscore=F)

overallTest(zscore50_5boot5.cull, metric.ls, tails =2,zscore=T)

t.tester(zscore50_5boot5.cull,metric.ls)

#genera where a difference in a network metric was observed
unique(zscore50_5.df$Genus[abs(zscore50_5.df$degree)>1.96])
unique(zscore50_5.df$Genus[abs(zscore50_5.df$d)>1.96])
unique(zscore50_5.df$Genus[abs(zscore50_5.df$weighted.closeness)>1.96])


######-------------------------------

#save and upload
save(zscore50_5.df,file="data/zscore50_5.RData")
#save(sexDiffsProp50_5.df,file='data/sexDiffsProp50_5.Rdata')

pb_upload("data/zscore50_5bootST.RData",
          name="zscore50_5bootST.RData",
          tag="data.v.1")




#bootstuff
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
