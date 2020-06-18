#3: Calculations and testing

#Runs analyses on the species-level network parameters calculated in
#2_netBuilding.R to compare males and females within the observed and
#simulated networks. 

################################################

#Species level network traits

################################################

###------------------
## Setup
pb_download("sex_trts_mix5.RData",
            dest="data",
            tag="data.v.1")

load("data/sex_trts_mix5.RData")

##specify the metrics I'll be looking at, number of cores to use
metric.ls <- c("degree","species.strength","weighted.betweenness",
               "weighted.closeness","d")
cores <- 10

##clean up the dataset, add some necessary columns to speed things up
traits5.ls <-
  mclapply(sex_trts_mix5, function(x){
  x$SiteYr <- paste(x$Site, x$Year, sep="_")
  x$Sp <- gsub( "_.*$", "", x$GenusSpecies )
  x <- filter(x,x$sex == "m" | x$sex == "f")
  droplevels(x$sex) #doesn't seem to work for whatever reason...
  return(x)
},mc.cores=cores)

###------------------
## calculate how different males and females are within each species
## within each SiteYr in each iteration. Comparison is 
## [male value - female value]: in the sexDiffs output, large values
## indicate that males had larger values than females.

sexDiffs5.df <- makeComp(traits5.ls, metric.ls, comparison = "diff")

## Saving and uploading after that long step
save(sexDiffs5.df, file = 'data/sexDiffs5.RData')

pb_upload('data/sexDiffs5.RData',
          name='sexDiffs5.RData',
          tag="data.v.1")

## **************************************************************** 
##Calculate how different the observed values were from the simulated values.
## when zscore=F, outupts the proportion of simulations where the value
## was less than or 50% equal to the obesrved difference value. 
## High proportions indicate that males were larger than females to a
## greater degree than the null expectation. Low proportions indicate that 
## females had higher values than males to a greater degree than expected

## downloading the output to skip above step 
pb_download("sexDiffs5.RData",
            dest="data",
            tag="data.v.1")
load('data/sexDiffs5.RData')

sexDiffsProp50_5.df <- calcNullProp50(sexDiffs5.df,
                                      metric.ls,
                                      zscore=FALSE)

##re-running above but to generate z-scores for graphs
zscore50_5.df <- calcNullProp50(sexDiffs5.df, 
                                metric.ls,
                                zscore=TRUE)
###------------------
##Test: proportion of species+sites where m v f difference in
##observed network was larger than many of the simulations. The output
##is the proportion of observations (Sp+Site+Year) where the male-female
## difference diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)

overallTest(sexDiffsProp50_5.df, metric.ls, tails=2, zscore=F)

######-------------------------------

#save and upload
save(zscore50_5.df,file="data/zscore50_5.RData")
save(sexDiffsProp50_5.df,file='data/sexDiffsProp50_5.Rdata')

pb_upload("data/zscore50_5.RData",
          name="zscore50_5.RData",
          tag="data.v.1")

pb_upload("data/sexDiffsProp50_5.Rdata",
          name="sexDiffsProp50_5.Rdata",
          tag="data.v.1")

# - Alternative code that will be removed before pub

######-------------------------------

sexDiffs2.df <- makeComp(traits2.ls, metric.ls, comparison = "diff")
sexDiffs3.df <- makeComp(traits3.ls, metric.ls, comparison = "diff")
sexDiffs10.df <- makeComp(traits10.ls, metric.ls, comparison = "diff")
sexDiffs20.df <- makeComp(traits20.ls, metric.ls, comparison = "diff")

#Add back in the GenusSpecies and Family columns 
sexDiffsProp50_5.df$GenusSpecies <- gsub("_.*$", "", sexDiffsProp50_5.df$SpSiteYr)
sexDiffsProp50_5.df$Family <- spec_all$Family[match(sexDiffsProp50_5.df$GenusSpecies,
                                                    spec_all$GenusSpecies)]

zscore50_5.df %>%
  mutate(GenusSpecies = gsub( "_.*$", "", SpSiteYr)) ->
  zscore50_5.df
zscore50_5.df$Family <- spec_all$Family[match(zscore50_5.df$GenusSpecies,
                                              spec_all$GenusSpecies)]
#add in a couple of missing families (NOTE: May not be necessary after data revision)
zscore50_5.df[432,8] <- "Syrphidae"
