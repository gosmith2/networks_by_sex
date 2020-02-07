## 1: Innitiation and dataset compiling

## Loads necessary libraries for the whole analysis process, downloads and cleans the original data from the Yosemite,
#Hedgerow, and Sky Islands, removing unused samples (e.g., those not caught with nets) and adds necessary columns. 
#It then compiles these datasets into a single object.

## setwd("~/Dropbox/networks_by_sex")

rm(list=ls())
library(piggyback) #
library(bipartite) #
library(tidyverse) #
library(parallel) #
library(nlme)

source('prepNets.R')
source('misc.R')

#Sys.getenv("GITHUB_PAT")
pb_download("specimens-yos.csv",
            dest="data",
            tag="data.v.1")

pb_download("specimens-hr.RData",
            dest="data",
            tag="data.v.1")

pb_download("specimens-si.RData",
            dest="data",
            tag="data.v.1")

####--------------------------####
#### YOSEMITE
####--------------------------####

spec.y <-read.csv("data/specimens-yos.csv")

## drop pan data
spec.y <- spec.y[spec.y$NetPan == "net",]

## drop extra round at L21 when field crew did not sample correctly
# from specimens
extra.round <- spec.y$Site == 'L21' & spec.y$Date == '2014-07-01'
spec.y <- spec.y[!extra.round,]

## get specimen data ready.
spec.y <- dat.clean(spec.y)
spec.y <- dat.dates(spec.y)

## drop non-bee, non-Syrphids
spec.y <- spec.y[spec.y$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
spec.y <- dat.rm.blanks(spec.y)
spec.y$GenusSpeciesSex<-ifelse(spec.y$Sex %in% c("m","f"),
                               paste(spec.y$GenusSpecies,
                                     spec.y$Sex,sep="_"),
                               paste(spec.y$GenusSpecies,
                                     "e",sep="_")
)



####--------------------------####
#### Hedgerow
####--------------------------####

load("data/specimens-hr.RData",verbose=TRUE)
#For whatever reason, resulting df is called "dd"
spec.h <- dd


## get specimen data ready.
spec.h <- dat.clean(spec.h)
spec.h$Date <- as.Date.character(spec.h$Date)
spec.h <- dat.dates(spec.h)


## drop non-bee, non-Syrphids
spec.h <- spec.h[spec.h$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
spec.h <- dat.rm.blanks(spec.h)

spec.h$GenusSpeciesSex <- ifelse(spec.h$Sex %in% c("m","f"),
                                 paste(spec.h$GenusSpecies,
                                       spec.h$Sex,sep="_"),
                                 paste(spec.h$GenusSpecies,
                                       "e",sep="_")
)

####--------------------------####
#### Sky Islands
####--------------------------####

load("data/specimens-si.RData",verbose=TRUE)

spec.s <- spec

spec.h <- spec.h[spec.h$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]
spec.s <- dat.rm.blanks(spec.s)

spec.s$GenusSpeciesSex <- ifelse(spec.s$Sex %in% c("m","f"),
                                 paste(spec.s$GenusSpecies,
                                       spec.s$Sex,sep="_"),
                                 paste(spec.s$GenusSpecies,
                                       "e",sep="_")
)

####--------------------------####
#### Combining network lists
####--------------------------####

keeps <- c("UniqueID",
           "GenusSpecies",
           "Site",
           "Year",
           "PlantGenusSpecies",
           "GenusSpeciesSex",
           "Sex")

bind_rows(select(spec.y,keeps),
          select(spec.h,keeps),
          select(spec.s,keeps)) -> spec.all

spec.all$SiteYr <- paste(spec.all$Site,spec.all$Year)

#save compiled dataset, upload to release
save(spec.all,file='data/spec.all.RData')

pb_upload("data/spec.all.RData",
          name="spec.all.RData",
          tag="data.v.1")