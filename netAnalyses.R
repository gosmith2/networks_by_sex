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

#trait.test<-traits.ls[1:3]

#bramlett 2014 and Bray2 2014 both have some NAs, even in [[1]]


makeLogRatio <- function(data, metrics, comparison="log") {
  ## takes the sex-level network traits of the pollinators and calculates
  ## the log ratio between males and females within each species and 
  ## network. Within each network, it also drops any species where only
  ## one sex was caught.
  
  som <- mclapply(data, function(x){
    
    #narrow down to each unique Site+Year combo
    col.ls <- lapply(unique(x$SiteYr),function(y){
      sites <- filter(x, x$SiteYr == y)
      
      #narrow down to each unique species within each SiteYr
      spp <- lapply(unique(sites$Sp), function(z){
        sp <- filter(sites,sites$Sp == z)
        
        #If there are two sexes of that species, get the log ratios of each 
        #metric specified in the metrics argument vector     
        if (length(sp$sex)==2) {
          ratios <- lapply(metrics, function(a){
              if(comparison = "log"){
                ratio <- log10(sp[,a][sp$sex == "m"]
                               /sp[,a][sp$sex == "f"])
              } else {
                diff <- sp[,a][sp$sex == "m"] - sp[,a][sp$sex == "f"]
              }
          #  }
            })
          #browser()
          ratios <- data.frame(ratios)
          colnames(ratios) <- metrics
          ratios$species <- z
          ratios$SiteYr <- y
        } 
        else {
            ratios <- NULL
        }
        return(ratios)
        #browser()
        })
      spp <- spp[!sapply(spp, is.null)]
      return(do.call(rbind, spp))
      browser()
      })
    #browser()
    return(do.call(rbind, col.ls))
  },mc.cores=cores)
  return(som)
}


#calculate how different males and females are in each iteration
logRatios.df <- makeLogRatio(traits.ls, metric.ls, comparison = "log")

sexDiffs.df <- makeLogRatio(traits.ls, metric.ls, comparison = "diff")

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

load("data/logRatios.RData")
load('data/sexDiffs.RData')

calcNullProp <- function(data, metrics, zscore=TRUE) {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    data[[x]]$SpSiteYr <- paste(data[[x]]$species,
                                data[[x]]$SiteYr,
                                sep="_")
    return(data[[x]])
  },mc.cores=cores)

  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)

  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(dist.df$SpSiteYr), function(y) {
    #browser()
    sp <- filter(dist.df, dist.df$SpSiteYr == y)
    obs <- filter(sp, sp$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      #browser()
      if (zscore == TRUE){
        metZ <- scale(sp[,z],center = TRUE, scale = TRUE)
        metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1]) 
        #gotta think about NA treatment here too. I kinda think its
        #getting rid of stuff again. Though this may be fixed when
        #I fix stuff above
      }else{
        metprop <- sum(sp[,z] <= obs[,z]) / length(sp$species)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SpSiteYr <- y
    return(mets)
  },mc.cores=cores)

  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}

zscores.ls <- calcNullProp(logRatios.df, metric.ls ,zscore=TRUE)

zscores.test<-calcNullProp(logRattest.df, metric.ls ,zscore=TRUE)

#these were long bits, so save above output for later
save(zscores.ls,file='data/zscores.RData')

pb_upload("data/zscores.RData",
          name="zscores.RData",
          tag="data.v.1")

pb_download("zscores.RData",
            dest="data",
            tag="data.v.1")

load("data/zscores.RData")



overallTest <- function(prop.dist, metrics, tails = 1, zscore = TRUE) {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2). Can alternatively
  ## give the proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time. 
  
  alpha <- lapply(metrics, function(x){
    #browser()
    clean <- !is.na(prop.dist[,x])
    clean.df <- prop.dist[clean,]
    if(tails == 1){
      #browser()
      if(zscore == TRUE) {
        sum(1.645 <= clean.df[,x]) / length(clean.df[,x])
      } else{
        sum(0.05 >= clean.df[,x]) / length(clean.df[,x])
      }
    } else {
      if(zscore == TRUE) {
        sum(1.96<=clean.df[,x] | -1.96>=clean.df[,x]) / length(clean.df[,x])
      } else{
        sum(0.025>=clean.df[,x] | 0.975<=clean.df) / length(clean.df[,x])
      }
    }
  })
  alpha <- data.frame(alpha)
  colnames(alpha) <- metrics
  return(alpha)
}

overallTest(zscores.ls,metric.ls,zscore=TRUE)
#returned: proportion of observations where m vs f differed. 
#near 0 = few differed significantly, near 1 = many differed significantly
  ## deg = 0.046, str = 0.056, weighted btw: NA, weighted close NA

spLevelTest <- function(prop.dist, metrics,zscore=TRUE, tails=1) {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2) for each
  ## species in the dataset separately. Can alternatively give the 
  ## proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time
  
  prop.dist$Sp <- gsub( "_.*$", "", prop.dist$SpSiteYr)
  spp <- lapply (unique(prop.dist$Sp),function(x){
    sp <- filter(prop.dist, prop.dist$Sp == x)
    sp.sig <- overallTest(sp, metrics, zscore=zscore, tails=tails)
    sp.sig$Sp <- x
    return(sp.sig)
    #browser()
  })
  spp.sig <- do.call(rbind,spp)
  #browser()
  return(spp.sig)
}

spLevelTest(zscores.ls,metric.ls)
#results: lots of zeros here too

#distribution of those z scores (1 value per sp per siteyr)
plot(density(zscores.ls$degree))
abline(v=1.645)
abline(v=c(-1.96,1.96))
  #its pretty clear that a TON of scores at exactly 0 (465/1063). 
  #I think we gotta assess what that means, how that should change 
  #things

plot(density(zscores.ls$species.strength))
abline(v=1.645)
  #307/1063 for species strength. 

##so why? Lets look at the actual traits
traits.ls[[1]] %>%
  filter(SiteYr == "MullerM_2010") %>%
  select(degree, Sp, sex)
  #Looks like many of the species that have both a male and a female
  #don't actually vary in their degree. For example, melisodes agilis
  #and lupina both had degrees of exactly 1 for both males and females
  #Presumably this is a consequence of small networks and/or many 
  #of the species specializing


#denisty plot with actual z distribution

genNullDist <- function(data, metrics, mean.by,zscore=TRUE) {
 
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #extract an average value w/in each siteyr across sims
  dist.build <- mclapply(unique(dist.df[,mean.by]), function(y) {
    #browser()
    net <- filter(dist.df, dist.df[,mean.by] == y)
    obs <- filter(net, net$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      #browser()
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

#so the above function works (right approach), but the distribution
#is essentially just around zero. like +/-0.3*10^-16. So the actual
#calculation is incorrect

nullDist.df <- genNullDist(logRatios.df,metric.ls,"sim",zscore=F)


#nullDist.test <- genNullDist(logRattest.df,metric.ls,"sim",zscore=F)

#this plot is right, i think. doesn't show z scores but looks good
plot(density(nullDist.df$degree,na.rm = T))
abline(v=-0.064)

meanDiffObs<-lapply(metric.ls, function(x){
  mean(logRatios.df[[1]][,x])
})


plotweb(nets.mix.clean[[1]]$MullerM.2010)


#regress amt of difference against absolute specialization? 
  #i.e, are more generalized sp more different b/w males and females?

#plotting: density shaded plot, single line at obs. just like she did on that paper

#for final project: above, density plot, make map w/ lat longs 



#######################################################

## Network level traits

#######################################################

pb_download("netlvlYH.RData",
            dest="data",
            tag="data.v.1")

load("data/netlvlYH.RData")

nullDistNets.df <- genNullDist(netlvl,metric.ls,"sim",zscore=F)











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




