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



####--------------------------####
#### YOSEMITE
####--------------------------####

#Sys.getenv("GITHUB_PAT")
#pb_download("specimens-yos.csv",
#            dest="data/specimens-yos.csv",
#            tag="data.v.1")
spec.y <-read.csv("data/specimens-yos.csv")

## drop pan data
spec.y <- spec.y[spec.y$NetPan == "net",]

## drop extra round at L21 when field crew did not sample correctly
# from specimens
extra.round <- spec.y$Site == 'L21' & spec.y$Date == '2014-07-01'
spec.y <- spec.y[!extra.round,]



## get specimen data ready. ##intuitive
dat.clean(spec.y)
dat.dates(spec.y)


## drop non-bee, non-Syrphids
spec.y <- spec.y[spec.y$Family %in% c("Andrenidae", "Apidae",
                                 "Colletidae", "Halictidae",
                                 "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
dat.rm.blanks(spec.y)
spec.y$GenusSpeciesSex<-paste(spec.y$GenusSpecies,spec.y$Sex)
spec$YearSR <- paste(spec$Year, spec$SampleRound, sep=".")

### site-level networks
build.nets(spec.y,"yos")
 
save(yos_GSSY,yos_GSY,yos_SSY,yos_SY, file="data/nets.Rdata")
#ah. saves the built networks elsewhere for analysis, which 
#will be similar to the analysis above



sp.lev <- calcSpec(nets, spec)

save(spec, file="../data/spec.Rdata")
write.csv(spec, file="../data/spec.csv", row.names=FALSE)


####--------------------------####
#### Hedgerow
####--------------------------####

#Sys.getenv("GITHUB_PAT")
#pb_download("specimens-yos.csv",
#            dest="data/specimens-yos.csv",
#            tag="data.v.1")
load("data/specimens_hr.RData",verbose=TRUE) 
  #For whatever reason, resulting df is called "dd"
spec.h<-dd


## get specimen data ready. ##intuitive
dat.clean(spec.h)
dat.dates(spec.h)

spec.h$PlantGenusSpecies <-  fix.white.space(paste(spec.h$PlantGenus,
                                                     spec.h$PlantSpecies,
                                                     spec.h$PlantVar,
                                                     spec.h$PlantSubSpecies))
spec.h$Int <-  fix.white.space(paste(spec.h$GenusSpecies,
                                       spec.h$PlantGenusSpecies))
spec.h$IntGen <-  fix.white.space(paste(spec.h$Genus,
                                          spec.h$PlantGenus))


## drop non-bee, non-Syrphids
spec.h <- spec.h[spec.h$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
dat.rm.blanks(spec.h)
spec.h$GenusSpeciesSex<-paste(spec.h$GenusSpecies,spec.h$Sex)
spec.h$YearSR <- paste(spec.h$Year, spec.h$SampleRound, sep=".")

### site-level networks
build.nets(spec.h,"hr")

save(hr_GSSY,hr_GSY,hr_SSY,hr_SY, file="data/nets.Rdata")



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


save(sp.lev, file='../data/splev.Rdata')


## checks
table(spec$GenusSpecies)

tab <- table(spec$GenusSpecies, spec$Site)
tab2 <- table(spec$PlantGenusSpecies, spec$Site)

table(spec$PlantGenusSpecies)
table(spec$PlantGenusSpecies, spec$Site)
table(spec$PlantGenusSpecies, spec$Year)

table(spec$PlantGenusSpecies, spec$Family)