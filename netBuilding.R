## Network (matrix) building from individual interactions. 
## Modified (very slightly) from Ponisio sky islands data prep.


##this should be run on lauren's computer via terminal (git bash): 
ssh gsmith@osmia.dyn.ucr.edu 
#ip is: 138.23.14.130
cd Documents
#git clone https://github.com/gosmith2/networks_by_sex.git
cd Documents/networks_by_sex

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

Sys.getenv("GITHUB_PAT")
pb_download("specimens-yos.csv",
            dest="specimens-yos.csv",
            tag="data.v.1")
pb_download("specimens-hr.RData",
            dest="specimens-hr.RData",
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
ssy.ls <- c(yos_SSY,hr_SSY)

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

cores <- 3

#randomize the sexes w/in species w/in sites
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
save(nets.mix.clean, file = 'data/mix_nets.RData')


#calculate network stats at the individual level, output into usable data frame
sex.trts.mix<-mclapply(nets.mix.clean,function(x) calcSpec(x), mc.cores = cores)

## confirm that the values are different
#sex.trts.mix[[1]] %>%
#	filter(Site == "Zamora", Year == 2014)

#write.csv(sex.trts.mix, file='data/sex_trts_mix.csv') #result here is 2.3 GB...
save(sex.trts.mix,file='data/sex_trts_mix.RData') #much better at 432 MB


#Sys.getenv("GITHUB_PAT")
pb_upload("data/sex_trts_mix.RData",
			name="sex_trts_mix.RData",
            tag="data.v.1")
pb_upload("data/mix_nets.RData",
			name="mix_nets.RData",
            tag="data.v.1")






###############################################

#below here is misc code that I don't run


##############################################


###############################################

#spec.all %>%
#  filter(SiteYr=="Zamora 2014")->
#  zam14.all

#ran.zam<-ran.sex(zam14.all)

#ram.zam <- do.call(rbind, unlist(ran.zam, recursive=FALSE))

#ram.zam %>%
#  filter(SiteYr=="Zamora 2014") %>%
#  (length(unique(sp$Sex)) != 1)
#  head()


sex_trts.df %>%
  filter(sex=="_") %>%
  select(GenusSpecies,Site)

#calculating network traits for males and females
sex_traits <- lapply(ssy.ls,function(x){
  y<-try(specieslevel(x))
  if(inherits(y, "try-error")) browser()
  return(y)
}
)


sex_trts.df
sex_trts.df$degree
all(sapply(ssy.ls,dim)>1) 


#remove the lower level (plants) to simplify indexing
pol_traits <- lapply(sex_traits,function(x){
  x[1]
    }
    )

poll_nodeSpecF <- lapply(sex_traits,function(x){
  lapply(x$'higher level',function(y){
    lapply(rownames(y),function(z){
      str_extract(z,"_.")
    }
    )
  }
  )
}
)
poll_nodeSpecF
    
poll_f <- ifelse(
  lapply(sex_traits,function(x){
    lapply(x$'higher level',function(y){
      lapply(rownames(y),function(z){
        str_extract(z,"_.")
    })})})=="_f",
  
    )
  }
  )
}
)

spec.h %>%
  filter(Sex!="f",Sex!="m") %>%
  select(Sex,GenusSpecies)

spec.y %>%
  filter(Genus=="Triepeolus") %>%
  select(Sex,GenusSpecies,GenusSpeciesSex)




#traceback. 
#sapply(nets, dim)<-give me all dimensions to see if any are actually vectors
#encircle area thats breaking with "try()"
  #then, "if(inherits(sl, "try-error")) browser()"



## *******************************************************************
## create a giant network to calculate specialization etc. acorss all
## SI
## *******************************************************************


#tst<-data.frame(SiteYr=c("a","a","a","b","b","c","c","c"),
#                ID=c(1:8),
#                GenusSpecies=c(1,1,2,1,1,3,3,3),
#                Sex=c("f","m","f","m","m","f","m","f"))




#tst$mix<-unlist(ran.sex(tst))

#tst.ls<-numeric()
#tst.ls[1]<-tst


#tst.ls<-ran.gen(tst,3)


#x$'higher level'$'node.specialisation.index.NSI'
#})


#tst<-sapply(rownames(sex_traits$Zamora.2014$'higher level'),function(x){
#  strsplit(x, split="_")
#})

#tst3<-sapply(rownames(sex_traits$Zamora.2014$'higher level'),function(x){
#  str_extract(x,"_.")
#})
#tst3

agg.spec <- aggregate(list(abund=spec$GenusSpecies),
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

#Drop Mariani 2012, MC1.2014, and H16.2012: too few interactions to calculate network parameters 
#m12<- spec.h$Site == "Mariani" & spec.h$Year == "2012"
#mc14<- spec.h$Site == "MC1" & spec.h$Year == '2014'
#h12<- spec.h$Site == "H16" & spec.h$Year == '2012'
#spec.h<-spec.h[!m12,]
#spec.h<-spec.h[!mc14,]
#spec.h<-spec.h[!h12,]




## checks
table(spec$GenusSpecies)

tab <- table(spec$GenusSpecies, spec$Site)
tab2 <- table(spec$PlantGenusSpecies, spec$Site)

table(spec$PlantGenusSpecies)
table(spec$PlantGenusSpecies, spec$Site)
table(spec$PlantGenusSpecies, spec$Year)

table(spec$PlantGenusSpecies, spec$Family)