## Network (matrix) building from individual interactions. 
## Modified (very slightly) from Ponisio sky islands data prep.


rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
library(tidyverse)
library(stringr)
library(nlme)
library(sciplot)
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
load("data/specimens_hr.RData",verbose=TRUE) 
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

#Drop Mariani 2012, MC1.2014, and H16.2012: too few interactions to calculate network parameters 
#m12<- spec.h$Site == "Mariani" & spec.h$Year == "2012"
#mc14<- spec.h$Site == "MC1" & spec.h$Year == '2014'
#h12<- spec.h$Site == "H16" & spec.h$Year == '2012'
#spec.h<-spec.h[!m12,]
#spec.h<-spec.h[!mc14,]
#spec.h<-spec.h[!h12,]


### site-level networks
build.nets(spec.h,"hr")

#save(hr_GSSY,hr_GSY,hr_SSY,hr_SY, file="data/hr_nets.Rdata")



####--------------------------####
#### Combining network lists
####--------------------------####

ssy.ls <- c(yos_SSY,hr_SSY)

#remove all networks with too few interactions to calculate metrics
ssy.ls <- ssy.ls[sapply(ssy.ls, function(x) all(dim(x) > 1))]

#calculate networks, output into usable data frame
sex_trts.df<-calcSpec(ssy.ls)


bargraph.CI(response=sex_trts.df$species.specificity.index,x.factor=sex_trts.df$sex)

nlme

sex_trts.df %>%
  filter(sex=="_") %>%
  select(GenusSpecies,Site)

spec.h %>%
  filter(Sex!="f",Sex!="m") %>%
  select(Sex,GenusSpecies)

spec.y %>%
  filter(Genus=="Triepeolus") %>%
  select(Sex,GenusSpecies,GenusSpeciesSex)




###############################################



##############################################


###############################################
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

  x$'higher level'$'node.specialisation.index.NSI'
})


tst<-sapply(rownames(sex_traits$Zamora.2014$'higher level'),function(x){
  strsplit(x, split="_")
})

tst3<-sapply(rownames(sex_traits$Zamora.2014$'higher level'),function(x){
  str_extract(x,"_.")
})
tst3


#traceback. 
#sapply(nets, dim)<-give me all dimensions to see if any are actually vectors
#encircle area thats breaking with "try()"
  #then, "if(inherits(sl, "try-error")) browser()"



## *******************************************************************
## create a giant network to calculate specialization etc. acorss all
## SI
## *******************************************************************
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


## checks
table(spec$GenusSpecies)

tab <- table(spec$GenusSpecies, spec$Site)
tab2 <- table(spec$PlantGenusSpecies, spec$Site)

table(spec$PlantGenusSpecies)
table(spec$PlantGenusSpecies, spec$Site)
table(spec$PlantGenusSpecies, spec$Year)

table(spec$PlantGenusSpecies, spec$Family)