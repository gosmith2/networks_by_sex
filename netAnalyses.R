#Runs analyses on the species-level network parameters calculated in
#netBuilding.R

rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
library(tidyverse)
library(piggyback)
library(stringr)
library(nlme)
library(sciplot)
library(parallel)
source('misc.R')
source('prepNets.R')

Sys.getenv("GITHUB_PAT")


################################################

#Species level network traits

################################################

###------------------
## Setup
pb_download("sex_trts_mixYH.RData",
            dest="data",
            tag="data.v.1")
pb_download("mix_netsYH.RData",
            dest="data",
            tag="data.v.1")

load("data/sex_trts_mixYH.RData")
load("data/mix_netsYH.RData") #object: nets.mix.clean

#specify the metrics I'll be looking at, number of cores to use
metric.ls <- c("degree","species.strength","weighted.betweenness","weighted.closeness" )

cores <- 10

##clean up the datasets, add some necessary columns to speed things up
traits.ls <- 
  mclapply(sex.trts.mix, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  droplevels(x$sex) #doesn't seem to work for whatever reason...
  return(x)
},mc.cores=cores)


###------------------
##calculate how different males and females within each species
##within each SiteYr are in each iteration

logRatios.df <- makeComp(traits.ls, metric.ls, comparison = "log")
  #log10 (male / female)

sexDiffs.df <- makeComp(traits.ls, metric.ls, comparison = "diff")
  #males - females


## Saving and uploading after that long step
save(logRatios.df, file = 'data/logRatios.RData')
save(sexDiffs.df, file = 'data/sexDiffs.RData')

pb_upload("data/logRatios.RData",
          name="logRatios.RData",
          tag="data.v.1")
pb_download("logRatios.RData",
            dest="data",
            tag="data.v.1")

pb_upload('data/sexDiffs.RData',
          name='sexDiffs.RData',
          tag="data.v.1")
pb_download("sexDiffs.RData",
            dest="data",
            tag="data.v.1")

load("data/logRatiosYH.RData")
load('data/sexDiffs.RData') #object: sexDiffs.df


###------------------
##Calculate how different the observed values were from the simulated

#by z-score
zscores.ls <- calcNullProp(logRatios.df, metric.ls ,zscore=TRUE)

#by proportion >= obs
sexDiffsProp.df <- calcNullProp(sexDiffs.df, metric.ls ,zscore=F) 

#by proportion > or 50% = obs
sexDiffsProp50.df <- calcNullProp50(sexDiffs.df,metric.ls, zscore=F)


## Saving and uploading after that long step
pb_upload("data/zscores.RData",
          name="zscores.RData",
          tag="data.v.1")

save(sexDiffsProp.df,file='data/sexDiffsPropYH.df')
save(sexDiffsProp50.df,file='data/sexDiffsProp50YH.df')


pb_upload("data/sexDiffsProp50YH.df",
          name="sexDiffsProp50YH.df",
          tag="data.v.1")

pb_download("zscores.RData",
            dest="data",
            tag="data.v.1")

pb_download("sexDiffsPropYH.df",
            dest="data",
            tag="data.v.1")

pb_download("sexDiffsProp50YH.df",
            dest="data",
            tag="data.v.1")

load("data/zscores.RData")
load("data/sexDiffsPropYH.df")
load("data/sexDiffsProp50YH.df")



###------------------
##Test: proportion of species+sites where m v f difference in 
##observed network was larger than many of the simulations 

overallTest(zscores.ls,metric.ls,zscore=TRUE)
#near 0 = few differed significantly, near 1 = many differed significantly
  ## deg = 0.046, str = 0.056, weighted btw: NA, weighted close NA

overallTest(sexDiffsProp.df,metric.ls,zscore=F)
  #

#same test as above, but with 50% >=, rather than >=
overallTest(sexDiffsProp50.df,metric.ls,zscore=F)


spLevelTest(sexDiffsProp50.df,metric.ls)
#results: lots of zeros here too


###--------------
## Generate null distributions and plotting 

#this one averages everything by simulation
nullDist.df <- genNullDist(logRatios.df,metric.ls,"sim",zscore=F)
save(nullDist.df,file='data/nullDist.RData')

#averages by sp+site+year (how overallTest was run)
nullDistDiff.df <- genNullDist(sexDiffs.df,metric.ls,zscore=F)

save(nullDistDiff.df,file='data/nullDistDiff.RData')


pb_upload("data/nullDist.RData",
          name="nullDist.RData",
          tag="data.v.1")

pb_download("nullDist.RData",
            dest="data",
            tag="data.v.1")

pb_upload("data/nullDistDiff.RData",
          name="nullDistDiff.RData",
          tag="data.v.1")

pb_download("nullDistDiff.RData",
            dest="data",
            tag="data.v.1")

#nullDist.test <- genNullDist(logRattest.df,metric.ls,"sim",zscore=F)

#meanObsZ<-lapply(metric.ls, function(x){
#  mean(logRatios.df[[1]][,x])
#})
#meanObsDiff<-lapply(metric.ls, function(x){
#  mean(sexDiffs.df[[1]][,x],na.rm=T)
#})

##this plot is right, i think. doesn't show z scores but looks good
#plot(density(nullDist.df$degree,na.rm = T))
#abline(v=-0.064)

##density plot for degree, absolute differences
plot(density(nullDistDiff.df$degree,na.rm = T))
abline(v=0)

plot(density(nullDistDiff.df$species.strength,na.rm = T))
abline(v=0)

plot(density(nullDistDiff.df$weighted.betweenness,na.rm = T))
abline(v=0)

plot(density(nullDistDiff.df$weighted.closeness,na.rm = T))
abline(v=0)


#regress amt of difference against absolute specialization? 
  #i.e, are more generalized sp more different b/w males and females?


#for final project: above, density plot, make map w/ lat longs 






#######################################################

## Network level traits

#######################################################

pb_download("netlvlYH.RData",
            dest="data",
            tag="data.v.1")

load("data/netlvlYH.RData")

overallTest(netlvl,metric.net,zscore=F)


metric.net <- c('connectance', 
                'number.of.compartments',
                'nestedness',
                'NODF',
                'robustness.HL',
                'robustness.LL',
                'vulnerability.LL',
                'H2',
                'niche.overlap.HL',
                'niche.overlap.LL',                  
                'functional.complementarity.HL',
                'functional.complementarity.LL')

calcNullPropNets <- function(data, metrics, zscore=TRUE) {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(dist.df$SiteYr), function(y) {
    #browser()
    site <- filter(dist.df, dist.df$SiteYr == y)
    obs <- filter(site, site$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      #browser()
      if (zscore == TRUE){
        metZ <- scale(site[,z],center = TRUE, scale = TRUE)
        #metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1]) 
        #gotta think about NA treatment here too. I kinda think its
        #getting rid of stuff again. Though this may be fixed when
        #I fix stuff above
      }else{
        metprop <- sum(site[,z] <= obs[,z]) / length(site$SiteYr)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SiteYr <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}

diffsNets.df <- calcNullPropNets(netlvl,metric.net,zscore=F) #quick

overallTest(diffsNets.df,metric.net,zscore=F)
  #yes, def some differences here. not all of them really really sig, but for sure
  


genNullDistNets <- function(data, metrics, mean.by,zscore=TRUE) {
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #extract an average value w/in each mean.by
  dist.build <- mclapply(unique(dist.df[,mean.by]), function(y) {
    net <- filter(dist.df, dist.df[,mean.by] == y)
    obs <- filter(net, net$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      browser()
      if (zscore == TRUE){
        mean.Z <- mean(scale(net[,z],center=T,scale=T),na.rm=T)
      }else{
        mean.value <- mean(net[,z])
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$mean.by <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,dist.build)
  return(sig.dist)
}

nullDistNets.df <- genNullDistNets(netlvl,metric.net,"sim",zscore=F)

save(nullDistNets.df,file="data/nullDistNetsYH.RData")

pb_upload("data/nullDistNetsYH.RData",
          name="nullDistNetsYH.RData",
          tag="data.v.1")
pb_download("nullDistNetsYH.RData",
          dest="data",
          tag="data.v.1")

meanObsDiffNet<-lapply(metric.net, function(x){
  mean(netlvl[[1]][,x],na.rm=T)
})

plot(density(nullDistNets.df$robustness.HL,na.rm = T))
abline(v=meanObsDiffNet[5])







####################################

##Below was exploratory, not run now

####################################


#bargraph.CI(response=sex_trts.df$d,x.factor=sex_trts.df$sex)

#Exploratory: 
##females higher degree (1.8 vs 1.45ish). same pattern diff #s for normalized
##**females higher str (0.4 vs. 0.25)
##males lower push pull (~-.7 vs -.6)
##males slightly higher nestedrank (0.55 vs 0.47ish)
##PDI about the same, v little var in either
##males slightly higher for resource.range (.92 vs .89? ish)
##sp. spec index about as above
##females higher PSI (0.19 vs .15)
##males v slightly higher node spec index (0.15ish vs 0.155?)
##**females much more between than males (0.055 vs 0.025). weighted even more diff
##females slightly more close than males (0.05 vs 0.047?). larger diff for weighted
##males higher fisher A (114 mil vs 96 mil)
##**females higher partner diversity (~0.3 vs 0.2)
##females higher eff. partners (1.58 vs 1.38)
##females slightly higher prop generality (0.4 vs 0.35)
##females slightly higher prop similarity (0.42 vs 0.38)
## d essentially the same




