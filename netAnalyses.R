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

makeLogRatio <- function(metrics) {
  som <- lapply(traits.ls, function(x){
    
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
tst.df<-makeLogRatio(metric)

tst.df[[1]] %>%
  filter(species == "Toxomerus marginatus") %>%
  filter(sex == "f" | sex == "m")

#take nulls
 # metrics, site, sex, sp. year, 
#pull out columns I want
#toss sp for which we don't have both males and females
#  log(m/f) for each site, yr, sp
}

calcNullProp <- function() {
  #make dist
  #calc proportion > or .5 equal to observed
  #outpus prop/1000 >equal site, yr, sex
  
  
}



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