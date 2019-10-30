## Network (matrix) building from individual interactions. 
## Modified (very slightly) from Ponisio sky islands data prep.


rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
source('misc.R')
source('prepNets.R')
#setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_prep.R')
#setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_make.R')
#setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/make_traditional.R')

#setwd("~/Dropbox/skyIslands/dataPrep")
#src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-read.csv("specimens-yos.csv")

#source("src/misc.R")
#source("src/prepNets.R")
#source("src/specialization.R")

## did not complete full sampling rounds in any of these sites. Was
## just scouting.
## can keep UK and SS when more species are IDed
## site.2.drop <- c("JM", "CC", "UK", "SS")
#site.2.drop <- c("JM", "CC", "SS")
#spec <- spec[!spec$Site %in% site.2.drop,]
#spec <- droplevels(spec)

## get specimen data ready. ##intuitive
spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                                           spec$Species,
                                           spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                                 spec$PlantSpecies,
                                                 spec$PlantVar,
                                                 spec$PlantSubSpecies))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies,
                                   spec$PlantGenusSpecies))
spec$IntGen <-  fix.white.space(paste(spec$Genus,
                                      spec$PlantGenus))

#set up dates. Intuitive
spec$Date <- as.Date(spec$Date, format='%m/%d/%y')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## drop non-bee, non-Syrphids
spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                 "Colletidae", "Halictidae",
                                 "Megachilidae", "Syrphidae"),]

##spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
##                                "Colletidae", "Halictidae",
##                                "Megachilidae"),]


## drop lasioglossum until we get 2017, 2018 IDs back from Joel
## spec <- spec[spec$Genus != "Lasioglossum",]


## for networks, drop specimens without plant IDs (or bee IDs)
spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]

## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds and there was basically nothing blooming the first
## round
#spec <- spec[!(spec$Site == "PL" & spec$Year == "2017"),]

## calculate orthoganol polynomials for doy. #??
spec$DoyPoly <- poly(spec$Doy, degree=2)
spec$DoyPoly1 <- spec$DoyPoly[,'1']
spec$DoyPoly2 <- spec$DoyPoly[,'2']
spec$DoyPoly <- NULL

## also for Latitude
spec$LatPoly <- poly(spec$Lat, degree=2)
spec$LatPoly1 <- spec$LatPoly[,'1']
spec$LatPoly2 <- spec$LatPoly[,'2']
spec$LatPoly <- NULL

save(spec, file="../data/spec.Rdata")
write.csv(spec, file="../data/spec.csv", row.names=FALSE)

## *******************************************************************
## create a giant network to calculate specialization etc. acorss all
## SI
## *******************************************************************
agg.spec <- aggregate(list(abund=spec$GenusSpecies), #not sure I understand the abund here
                      list(GenusSpecies=spec$GenusSpecies,
                           PlantGenusSpecies=spec$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                          agg.spec$GenusSpecies,
                          agg.spec$abund, FUN=sum)

all.traits <- specieslevel(nets.all)
## calculate rarified plant.pol degree
rare.plants.degree <- apply(nets.all, 1, chao1)
rare.pols.degree <- apply(nets.all, 2, chao1)

traits <- data.frame(GenusSpecies= unlist(sapply(all.traits,
                                                 rownames)),
                     do.call(rbind, all.traits))

traits$r.degree <-  rare.pols.degree[match(traits$GenusSpecies,
                                           names(rare.pols.degree))]
traits$r.degree[is.na(traits$r.degree)] <-
  rare.plants.degree[match(traits$GenusSpecies[is.na(traits$r.degree)],
                           names(rare.plants.degree))]

rownames(traits) <- NULL

write.csv(traits, file='../data/traits.csv')

### site-level networks
spec$YearSR <- paste(spec$Year, spec$SampleRound, sep=".")

nets <- breakNet(spec, 'Site', 'YearSR') #all of the networks, split by $ somehow. does all the aggregating etc
graphs <- lapply(nets, graph.incidence, weighted=TRUE) #makes a graph for each?
save(graphs, nets, file="../data/nets.Rdata")
  #ah. saves the built networks elsewhere for analysis, which 
  #will be similar to the analysis above
 


sp.lev <- calcSpec(nets, spec)
save(sp.lev, file='../data/splev.Rdata')


## checks
table(spec$GenusSpecies)

tab <- table(spec$GenusSpecies, spec$Site)
tab2 <- table(spec$PlantGenusSpecies, spec$Site)

table(spec$PlantGenusSpecies)
table(spec$PlantGenusSpecies, spec$Site)
table(spec$PlantGenusSpecies, spec$Year)

table(spec$PlantGenusSpecies, spec$Family)