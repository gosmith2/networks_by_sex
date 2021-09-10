## 1: Innitiation and dataset compiling

## Loads necessary libraries for the whole analysis process, downloads and cleans the original data from the Yosemite,
#Hedgerow, and Sky Islands datasets, removing unused samples (e.g., those not caught with nets) and adds necessary columns. 
#It then compiles these datasets into a single object.

## setwd("~/Dropbox/networks_by_sex")

rm(list=ls())
library(piggyback) 
library(bipartite)
library(tidyverse) 
library(parallel) 
library(vegan)
library(stringr)
library(nlme)
library(fossil)
library(reshape2)

source('prepNets.R')


##--------Authenticate so that you can download large datasets with piggyback:

#Copy the string below without the hash:
# GITHUB_PAT=7423a37de1b5d3913210273a64d0d89e65e5c8f9
  #This is the git token that will give piggyback access to the data releases.

#open your r environment document
usethis::edit_r_environ()

  #Paste the string you coppied with the git token. 
  #Press enter to leave 1 blank space below the token
  #save and close the environment doc. 
  #Restart R and run from below

Sys.getenv("GITHUB_PAT")

#For more info throughout this, see
#https://cran.r-project.org/web/packages/piggyback/vignettes/intro.html



#----Download datasets

pb_download("specimens-yos.csv",
            dest="data",
            tag="data.v.1")

pb_download("specimens-hr.RData",
            dest="data",
            tag="data.v.1")

pb_download("specimens-si.csv",
            dest="data",
            tag="data.v.1")

####--------------------------####
#### YOSEMITE
####--------------------------####

spec.y <-read.csv("data/specimens-yos.csv")

## drop pan data
spec.y <- spec.y[spec.y$NetPan == "net",]

## drop extra round at L21 when field crew did not sample correctly
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
## add a column combining sites and years
spec.y$SiteYr <- paste(spec.y$Site,spec.y$Year)

## add a column for dataset source
spec.y$dataset <- rep("y")

save(spec.y,file='data/spec_y.RData')


####--------------------------####
#### Hedgerow
####--------------------------####

load("data/specimens-hr.RData",verbose=TRUE)
#resulting df is called "dd"
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

## add a column combining sites and years
spec.h$SiteYr <- paste(spec.h$Site,spec.h$Year)

## add a column for dataset source
spec.h$dataset <- rep("h")

## remove observations for which the pollinator was only identified to Genus
spec.h <- spec.h[str_detect(spec.h$GenusSpecies," "),]

save(spec.h,file='data/spec_h.RData')

####--------------------------####
#### Sky Islands
####--------------------------####

spec.s <- read.csv("data/specimens-si.csv")

spec.s <- spec.s[spec.s$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]

## drop pan data
spec.s <- spec.s[spec.s$Method == "Net",]

#combine plant and pollinator names into single columns, do other cleanup
spec.s <- dat.clean(spec.s)


#remove blanks
spec.s <- dat.rm.blanks(spec.s)

#add some more necessary columns
spec.s$GenusSpeciesSex <- ifelse(spec.s$Sex %in% c("m","f"),
                                 paste(spec.s$GenusSpecies,
                                       spec.s$Sex,sep="_"),
                                 paste(spec.s$GenusSpecies,
                                       "e",sep="_")
)

spec.s$Year <- as.numeric(paste0("20",str_sub(spec.s$Date,-2)))

## add a column combining sites and years
spec.s$SiteYr <- paste(spec.s$Site,spec.s$Year)

## add a column for dataset source
spec.s$dataset <- rep("s")

save(spec.s,file='data/spec_s.RData')


####--------------------------####
#### Combining network lists
####--------------------------####

keeps <- c("UniqueID",
           "GenusSpecies",
           "Site",
           "Year",
           "PlantGenusSpecies",
           "GenusSpeciesSex",
           "Sex",
           "SiteYr",
           "Family",
           "dataset")

#all sites
bind_rows(select(spec.y,keeps),
          select(spec.h,keeps),
          select(spec.s,keeps)) -> spec_all

#save compiled dataset, upload to release
save(spec_all,file='data/spec_all.RData')

pb_upload("data/spec_all.RData",
          name="spec_all.RData",
          tag="data.v.1")
