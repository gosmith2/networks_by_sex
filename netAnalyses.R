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
source('misc.R')
source('prepNets.R')

#Sys.getenv("GITHUB_PAT")
#pb_download("sex_trts_mix.RData",
#            dest="data",
#            tag="data.v.1")
#pb_download("mix_nets.RData",
#            dest="data",
#            tag="data.v.1")

load("data/sex_trts_tst.RData")
load("data/mix_nets.RData")


#clean up the datasets, add some necessary columns to speed things up
traits.ls <- lapply(sex.trts.mix, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  return(x)
})
trait.test<-traits.ls[1:3]

makeLogRatio <- function(data, metrics) {
  som <- lapply(data, function(x){
    
    #narrow down to each unique Site+Year combo
    col.ls <- lapply(unique(x$SiteYr),function(y){
      sites <- filter(x, x$SiteYr == y)
      
      #narrow down to each unique species within each SiteYr
      spp <- lapply(unique(sites$Sp), function(z){
        sp <- filter(sites,sites$Sp == z)
        
        #If there are two sexes of that species, get the log ratios of each 
        #metric specified in the metrics argument vector     
        if (length(sp$sex)==2) {
          ratios<-lapply(metrics, function(a){
            ratio <- log10(sp[,a][sp$sex == "m"]/sp[,a][sp$sex == "f"])
          
            })
          #browser()
          ratios <- data.frame(ratios)
          colnames(ratios) <- metrics
          ratios$species <- z
          ratios$SiteYr <-y
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
  })
  return(som)
}
metric<-c("degree","species.strength")
tst.df<-makeLogRatio(traits.ls, metric)


calcNullProp <- function(data, metrics, zscore=TRUE) {
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- lapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    data[[x]]$SpSiteYr <- paste(data[[x]]$species,
                                data[[x]]$SiteYr,
                                sep="_")
    return(data[[x]])
  })
  #browser()
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)

  #combine values for sp+yr+site
  sigLevel <- lapply(unique(dist.df$SpSiteYr), function(y) {
    #browser()
    sp <- filter(dist.df, dist.df$SpSiteYr == y)
    obs <- filter(sp, sp$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      #browser()
      if (zscore == TRUE){
        metZ <- scale(sp[,z],center = TRUE, scale = TRUE)
        metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1])
      }else{
        metprop <- sum(sp[,z] <= obs[,z]) / length(sp$species)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SpSiteYr <- y
    return(mets)
  })
  #browser()
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}

tst.sig <- calcNullProp(tst.df,metric,zscore=TRUE)

overallTest <- function(prop.dist, metrics, tails = 1, zscore = TRUE) {
  alpha <- lapply(metrics, function(x){
    if(tails == 1){
      if(zscore == TRUE) {
        sum(1.645 <= prop.dist[,x]) / length(prop.dist[,x])
      } else{
        sum(0.05 >= prop.dist[,x]) / length(prop.dist[,x])
      }
    } else {
      if(zscore == TRUE) {
        sum(1.96<=prop.dist[,x] | -1.96>=prop.dist[,x]) / length(prop.dist[,x])
      } else{
        sum(0.025>=prop.dist[,x] | 0.975<=prop.dist) / length(prop.dist[,x])
      }
    }
  })
  alpha <- data.frame(alpha)
  colnames(alpha) <- metrics
  return(alpha)
}

overallTest(tst.sig,metric,zscore=TRUE)

#add a toggle in for z scores (which would be good as some sort of measure of effect size between metrics)
#plotting: density shaded plot, single line at obs. just like she did on that paper

#for final project: above, density plot, make map w/ lat longs 

spLevelTest <- function(prop.dist, metrics,zscore=TRUE, tails=1) {
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

spLevelTest(tst.sig,metric)



head(sex_trts.df)

bargraph.CI(response=sex_trts.df$d,x.factor=sex_trts.df$sex)

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

#closeness, betweenness, degree, str




#fxn "sample", loop (apply) over sites, over species, then reorganize sex based on the vector
#end goal: probably giant heirarcical list. each element is essentially equivalent to ssy.ls or sex_trts.df. 
  #should be 1000 total: 999 random + real. 


#mclapply: can run a whole bunch of things in parallel to run through the stuff quickly





nlme

sex_trts.df %>%
  filter(sex=="_") %>%
  select(GenusSpecies,Site)