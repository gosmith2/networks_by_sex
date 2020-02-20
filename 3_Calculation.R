#Runs analyses on the species-level network parameters calculated in
#netBuilding.R

################################################

#Species level network traits

################################################

###------------------
## Setup
pb_download("sex_trts_mixYHS2.RData",
            dest="data",
            tag="data.v.1")
pb_download("mix_netsYHS.RData",
            dest="data",
            tag="data.v.1")

load("data/sex_trts_mixYHS2.RData")
load("data/mix_netsYH.RData") #object: nets.mix.clean

#specify the metrics I'll be looking at, number of cores to use
metric.ls <- c("degree","species.strength","weighted.betweenness",
               "weighted.closeness","d")

cores <- 10

##clean up the datasets, add some necessary columns to speed things up
traits20.ls <-
  mclapply(sex.trts.mix20, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  droplevels(x$sex) #doesn't seem to work for whatever reason...
  return(x)
},mc.cores=cores)


###------------------
##calculate how different males and females within each species
##within each SiteYr are in each iteration

sexDiffs2.df <- makeComp(traits2.ls, metric.ls, comparison = "diff")

sexDiffs3.df <- makeComp(traits3.ls, metric.ls, comparison = "diff")
sexDiffs5.df <- makeComp(traits5.ls, metric.ls, comparison = "diff")
sexDiffs10.df <- makeComp(traits10.ls, metric.ls, comparison = "diff")
sexDiffs20.df <- makeComp(traits20.ls, metric.ls, comparison = "diff")


#17/196 networks are from SI, rest from HR. NONE from Yos


## Saving and uploading after that long step
save(sexDiffs20.df, file = 'data/sexDiffs20.RData')

pb_upload('data/sexDiffs20.RData',
          name='sexDiffs20.RData',
          tag="data.v.1")

pb_download("sexDiffs20.RData",
            dest="data",
            tag="data.v.1")

load('data/sexDiffs20.RData')

###------------------
##Calculate how different the observed values were from the simulated

sexDiffsProp50_20.df <- calcNullProp50(sexDiffs20.df,
                                      metric.ls,
                                      zscore=F)

zscore50_2.df <- calcNullProp50(sexDiffs2.df, metric.ls)

zscore50_2.df %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  zscore50_2.df

sexDiffsProp50_5.df$GenusSpecies <- gsub("_.*$", "", sexDiffsProp50_5.df$SpSiteYr)
sexDiffsProp50_5.df$Family <- spec.all$Family[match(sexDiffsProp50_5.df$GenusSpecies,
                                              spec.all$GenusSpecies)]

zscore50_2.df$Family <- spec.all$Family[match(zscore50_2.df$GenusSpecies,
                                              spec.all$GenusSpecies)]
zscore50_2.df[432,8] <- "Syrphidae"

save(zscore50_2.df,file="data/zscore50_2.RData")

save(sexDiffsProp50_20.df,file='data/sexDiffsProp50YH_20.Rdata')

pb_upload("data/zscore50_2.RData",
          name="zscore50_2.RData",
          tag="data.v.1")

pb_upload("data/sexDiffsProp50YH_20.Rdata",
          name="sexDiffsProp50YH_20.Rdata",
          tag="data.v.1")



