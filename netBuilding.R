## Network (matrix) building from individual interactions. 
## Modified (very slightly) from Ponisio sky islands data prep.


##this should be run on lauren's computer via terminal (git bash): 
#ssh gsmith@osmia.dyn.ucr.edu 
#ip is: 138.23.14.130
#cd Documents
#git clone https://github.com/gosmith2/networks_by_sex.git
#cd Documents/networks_by_sex

rm(list=ls())
library(piggyback)
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
library(tidyverse)
library(stringr)
library(nlme)
library(sciplot)
library(parallel)

source('misc.R')
source('prepNets.R')

#Sys.getenv("GITHUB_PAT")
pb_download("specimens-yos.csv",
            dest="data",
            tag="data.v.1")
pb_download("specimens-hr.RData",
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
spec.y<-dat.clean(spec.y)
spec.y<-dat.dates(spec.y)


## drop non-bee, non-Syrphids
spec.y <- spec.y[spec.y$Family %in% c("Andrenidae", "Apidae",
                                 "Colletidae", "Halictidae",
                                 "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
spec.y<-dat.rm.blanks(spec.y)
spec.y$GenusSpeciesSex<-ifelse(spec.y$Sex %in% c("m","f"),
                               paste(spec.y$GenusSpecies,
                                     spec.y$Sex,sep="_"),
                               paste(spec.y$GenusSpecies,
                                     "e",sep="_")
)
spec.y$YearSR <- paste(spec.y$Year, 
                       spec.y$SampleRound, 
                       sep=".")

### site-level networks
build.nets(spec.y,"yos") #uses the breakNets and breakNetsSex fxns to
#aggregate and then samp2site larger network into chunks
 
#save(yos_GSSY,yos_GSY,yos_SSY,yos_SY, file="data/yos_nets.Rdata")

#sp.lev <- calcSpec(nets, spec)

#save(spec, file="../data/spec.Rdata")
#write.csv(spec, file="../data/spec.csv", row.names=FALSE)


####--------------------------####
#### Hedgerow
####--------------------------####

#Sys.getenv("GITHUB_PAT")
#pb_download("specimens-hr.RData",
#            dest="data/specimens-hr.RData",
#            tag="data.v.1")
#may have to do the following on Lauren's computer (Osmia):
load("data/specimens-hr.RData",verbose=TRUE) 
  #For whatever reason, resulting df is called "dd"
spec.h<-dd


## get specimen data ready.
spec.h<-dat.clean(spec.h)
spec.h$Date<-as.Date.character(spec.h$Date)
spec.h<-dat.dates(spec.h)



## drop non-bee, non-Syrphids
spec.h <- spec.h[spec.h$Family %in% c("Andrenidae", "Apidae",
                                      "Colletidae", "Halictidae",
                                      "Megachilidae", "Syrphidae"),]

## for networks, drop specimens without plant IDs (or bee IDs)
spec.h<-dat.rm.blanks(spec.h)

spec.h$GenusSpeciesSex<-ifelse(spec.h$Sex %in% c("m","f"),
                               paste(spec.h$GenusSpecies,
                                     spec.h$Sex,sep="_"),
                               paste(spec.h$GenusSpecies,
                                     "e",sep="_")
)

spec.h$YearSR <- paste(spec.h$Year,
                       spec.h$SampleRound, 
                       sep=".")


### site-level networks
build.nets(spec.h,"hr")

#save(hr_GSSY,hr_GSY,hr_SSY,hr_SY, file="data/hr_nets.Rdata")



####--------------------------####
#### Combining network lists
####--------------------------####

#observed networks
#ssy.ls <- c(yos_SSY,hr_SSY)

##randomizing males and females for nulls
keeps<-c("UniqueID",
          "GenusSpecies",
          "Site",
          "Year",
          "PlantGenusSpecies",
          "GenusSpeciesSex",
          "Sex")

bind_rows(select(spec.y,keeps),
          select(spec.h,keeps))->
  spec.all

spec.all$SiteYr<-paste(spec.all$Site,spec.all$Year)

#cores should be 3 for bombus, 10 for osmia
cores <- 10

#randomize the sexes w/in species w/in sites. observed network is element 1
rand.sexes.ls<-ran.gen(spec.all,999,cores)

## checking that they mixed
#rand.sexes.ls[[1]] %>%
#	filter(SiteYr=='Zamora 2014') %>%
#	select(GenusSpeciesMix,UniqueID)

#build the networks at the sex level using the mixed sexes
nets.mix<-mclapply(rand.sexes.ls,function(y){
  breakNetMix(y,'Site','Year','GenusSpeciesMix')
}, mc.cores = cores)

## confirm that the networks are actually different
#nets.mix[[1]]$Zamora.2014
#nets.mix[[2]]$Zamora.2014

#remove all networks with too few interactions to calculate metrics
nets.mix.clean<-mclapply(nets.mix, function(x){
  x[sapply(x,function(y) all(dim(y)>1))]
}, mc.cores=cores)

#save the networks themselves
save(nets.mix.clean, file = 'data/mix_netsYH.RData')


#calculate network stats at the individual level, output into usable data frame
sex.trts.mix<-mclapply(nets.mix.clean,function(x) calcSpec(x), mc.cores = cores)

## confirm that the values are different
#sex.trts.mix[[1]] %>%
#	filter(Site == "Zamora", Year == 2014)

save(sex.trts.mix,file='data/sex_trts_mixYH.RData')


#Sys.getenv("GITHUB_PAT")
pb_upload("data/sex_trts_mixYH.RData",
			name="sex_trts_mixYH.RData",
            tag="data.v.1")
pb_upload("data/mix_netsYH.RData",
			name="mix_netsYH.RData",
            tag="data.v.1")






###############################################

#below here is misc code that I don't run


##############################################


###############################################



#traceback. 
#sapply(nets, dim)<-give me all dimensions to see if any are actually vectors
#encircle area thats breaking with "try()"
  #then, "if(inherits(sl, "try-error")) browser()"



## checks
table(spec$GenusSpecies)

tab <- table(spec$GenusSpecies, spec$Site)
tab2 <- table(spec$PlantGenusSpecies, spec$Site)

table(spec$PlantGenusSpecies)
table(spec$PlantGenusSpecies, spec$Site)
table(spec$PlantGenusSpecies, spec$Year)

table(spec$PlantGenusSpecies, spec$Family)