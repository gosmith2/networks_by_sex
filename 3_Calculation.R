#3: Calculations

#Runs analyses on the species-level network parameters calculated in
#2_netBuilding.R to compare males and females within the observed and
#simulated networks. 


library(parallel)
library(piggyback)
library(bipartite)
library(tidyverse)
source('prepNets.R')

#download and load the observed networks
pb_download("nets_mix_clean10kO.RData",
            dest="data",
            tag="data.v.1")
load("data/nets_mix_clean10kO.RData")

#download + load the validation networks
pb_download("nets_mix_clean10kS.RData",
            dest="data",
            tag="data.v.1")
load("data/nets_mix_clean10kS.RData")

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

#save
save(sex_trts_mix5O,file='data/sex_trts_mix5O.RData')
save(sex_trts_mix5S,file='data/sex_trts_mix5S.RData')


##specify the number of cores to use
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


#pb_upload('data/sexDiffs5S.RData',
#          name='sexDiffs5S.RData',
#          tag="data.v.1")

#pb_upload('data/sexDiffs5O.RData',
#          name='sexDiffs5O.RData',
#          tag="data.v.1")


## **************************************************************** 
##Calculate how different the observed values were from the simulated values.
## when zscore=F, outupts the proportion of simulations where the value
## was less than or 50% equal to the obesrved difference value. 
## High proportions indicate that males were larger than females to a
## greater degree than the null expectation. Low proportions indicate that 
## females had higher values than males to a greater degree than expected

## downloading the output to skip above calculation steps 
#pb_download("sexDiffs5O.RData",
#            dest="data",
#            tag="data.v.1")
#load('data/sexDiffs5O.RData')


##generate z-scores
zscore50_O.df <- calcNullProp50(sexDiffs5O.df, 
                                metric2.ls,
                                zscore=TRUE)
zscore50_S.df <- calcNullProp50(sexDiffs5S.df, 
                                   metric2.ls,
                                   zscore=TRUE)

save(zscore50_O.df,file='data/zscore50_O.RData')
save(zscore50_S.df,file='data/zscore50_S.RData')


#pb_upload('data/zscore50_S.RData',name='zscore50_S.RData',tag='data.v.1')
#pb_upload('data/zscore50_O.RData',name='zscore50_O.RData',tag='data.v.1')

