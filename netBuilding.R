## Network (matrix) building from individual interactions. 
## Modified (very slightly) from Ponisio sky islands data prep.


##this should be run on lauren's computer via terminal (git bash): 
#ssh gsmith@osmia.dyn.ucr.edu 
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
library(SYNCSA)

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


####--------------------------####
#### Hedgerow
####--------------------------####

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



####--------------------------####
#### Combining network lists
####--------------------------####

keeps<-c("UniqueID",
          "GenusSpecies",
          "Site",
          "Year",
          "PlantGenusSpecies",
          "GenusSpecies",
          "GenusSpeciesSex",
          "Sex")

bind_rows(select(spec.y,keeps),
          select(spec.h,keeps))->
  spec.all

spec.all$SiteYr<-paste(spec.all$Site,spec.all$Year)



####-------------------------------------------####
#### Randomizing m and f for specieslevel analyses
####-------------------------------------------####


#cores should be 3 for bombus, 10 for osmia
cores <- 10

#randomize the sexes w/in species w/in sites. observed network is element 1
rand.sexes.ls<-ran.gen(spec.all,999,cores)

## checking that they mixed
#rand.sexes.ls[[1]] %>%
#	filter(SiteYr=='Zamora 2014') %>%
#	select(GenusSpeciesMix,UniqueID)

#build the networks at the sex level using the mixed sexes
nets.mix<-mclapply(rand.sexes.ls, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
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

nets.mix.test <- nets.mix.clean[1:3]

#calculate network stats at the individual level, output into usable data frame
sex.trts.mix<-mclapply(nets.mix.clean,function(x) calcSpec(x), mc.cores = cores)

## confirm that the values are different
#sex.trts.mix[[1]] %>%
#	filter(Site == "Zamora", Year == 2014)

save(sex.trts.mix,file='data/sex_trts_mixYH.RData')



#calculate network-level traits (likely similar b/w obs and sim networks, but checking)

#netlvl <- mclapply(nets.mix.clean, function(x){
#  #browser()
#  network.lvl <- lapply(names(x), function(y){
#    #browser()
#    nl <- networklevel(x[[y]])
#    #nl <- data.frame(nl)
#    #nl$SiteYr <- y
#    })
#  
#  traits <- data.frame(do.call(rbind, network.lvl))
#  traits$SiteYr <- names(x)
#  return(traits)
#},mc.cores=cores)

#save(netlvl,file="data/netlvlYH.RData")


#Sys.getenv("GITHUB_PAT")
pb_upload("data/sex_trts_mixYH.RData",
			name="sex_trts_mixYH.RData",
            tag="data.v.1")
pb_upload("data/mix_netsYH.RData",
			name="mix_netsYH.RData",
            tag="data.v.1")
pb_upload("netlvlYH.RData",
          name="netlvlYH.RData",
          tag="data.v.1")



####----------------------------####
#### Prep for networklevel analyses
####----------------------------####

metric.net <- c('connectance', 
                'number of compartments',
                'nestedness',
                'NODF',
                'H2',
                'robustness',
                'vulnerability',
                'niche overlap',
                'functional complementarity')

N = 999
cores = 10

##all networks, with males and females included separately
nets.obs.sex<-breakNetMix(spec.all,'Site','Year',"GenusSpeciesSex")

netStats.sex <- mclapply(nets.obs.sex, 
                      calcNetworkMetrics, 
                      N=N,
                      index=metric.net,
                      mc.cores=cores)

save(netStats.sex, file = "data/netStats_sex.RData")

pb_upload("data/netStats_sex.RData",
          name="netStats_sex.RData",
          tag="data.v.1")

##all networks, males and females lumped into species
nets.obs.sp <- breakNet(spec.all,'Site','Year')
netStats.sp <- mclapply(nets.obs.sp, 
                         calcNetworkMetrics, 
                         N=N,
                         index=metric.net,
                         mc.cores=cores)

save(netStats.sp, file = "data/netStats_sp.RData")

pb_upload("data/netStats_sp.RData",
          name="netStats_sp.RData",
          tag="data.v.1")

##all networks, males dropped
spec.all.drop <- filter(spec.all, spec.all$Sex=="f")
nets.obs.f <- breakNet(spec.all.drop,'Site','Year')

netStats.f <- mclapply(nets.obs.f, 
                        calcNetworkMetrics, 
                        N=N,
                        index=metric.net,
                        mc.cores=cores)

save(netStats.f, file = "data/netStats_f.RData")

pb_upload("data/netStats_f.RData",
          name="netStats_f.RData",
          tag="data.v.1")


##Load and merge the 3 stats datasets
pb_download("netStats_f.RData",
          dest="data",
          tag="data.v.1")
pb_download("netStats_sex.RData",
            dest="data",
            tag="data.v.1")
pb_download("netStats_sp.RData",
            dest="data",
            tag="data.v.1")

load("data/netStats_sp.RData")
load("data/netStats_f.RData")
load("data/netStats_sex.RData")


names(netStats.f)

netStats.all <- do.call(rbind,
                        lapply(names(netStats.f), 
                               function (x) {
                                 
  all.df <- data.frame(rbind(netStats.sex[x][[1]][1:39],
                             netStats.sp[x][[1]][1:39],
                             netStats.f[x][[1]][1:39]))
  names(all.df) <- names(netStats.sex["Zamora.2014"][[1]])
  all.df$SiteYr <- x
  all.df$trt <- c("sex","sp","fem") 
  

  return(all.df)
})
)

ggplot(netStats.all, aes(x=zrobustness.HL, color=trt)) +
  geom_density()
