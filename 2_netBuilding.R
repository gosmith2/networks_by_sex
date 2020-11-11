## 2: Network (matrix) building from individual interactions.

## Loads the compiled dataset, generates simulated data randomizing sex,
## then breaks the observations into network matrices at the Site+Year level. 
## Lastly, it calculates node-level network metrics
## Modified (very slightly) from Ponisio sky islands data prep.

## all working directories are relative to the gitrepo 
library(parallel)
library(piggyback)
source("prepNets.R")

pb_download("spec_all.RData",
            dest="data",
            tag="data.v.1")

load("data/spec_all.RData")

####-------------------------------------------####
#### Randomizing males and females for analyses
####-------------------------------------------####

## WARNING: the following code is computationally intesnive.
## It is currently written to be processed by 10 cores
## via the function mclapply. 


## the number of cores you wish to run the simulation in parallel across 
cores <- 10

## randomize the sexes w/in species w/in sites. observed data is
## element 1 of the resulting list, randomized are elements 2-1000

rand_sexes2k <- ran.gen(spec_all, 1999, cores)
  #Note: 2000 was the highest number of iterations that successfully compiled, 
  #presumably due to memory constraints

save(rand_sexes2k,file="data/rand_sexes2k.RData")

#build the networks at the sex level using the mixed sexes
nets_mix <- mclapply(rand_sexes2k, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)


## confirm that the networks are actually different
ifelse(any(nets_mix[[1]][[1]]!=nets_mix[[2]][[1]]),
       print("SUCCESS: the networks are different"),
       print("WARNING: the networks are not different"))

#remove all networks with too few interactions to calculate metrics
nets_mix_clean2k <- mclapply(nets_mix, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

#save the networks themselves and upload them
save(nets_mix_clean2k, file = 'data/nets_mix_clean2k.Rdata')
pb_upload("data/nets_mix_clean2k.RData",
          name="nets_mix_clean2k.RData",
          )
pb_download('nets_mix_clean2k.RData',
            dest='data',
            tag="data.v.1")

####-------------------------------------------####
#### Calculate node-level network parameters
####-------------------------------------------####

## calculate network stats at the node level, output into usable data frame
# reshuffling threshold is 5 in this case (i.e. there must be at least 5
# males and 5 females of a given species in a given network for it to be included)
sex_trts_mix5 <- mclapply(nets_mix_clean2k,
                          function(x) calcSpec(x, indiv = 5),
                          mc.cores = cores)

## confirm that the values are different
ifelse(any(sex_trts_mix5[[1]]$weighted.closeness!=
             sex_trts_mix5[[2]]$weighted.closeness),
       print("SUCCESS: the network statistics are different between randomizations"),
       print("WARNING: the network statistics are not different between randomizations")
)

#save and upload
save(sex_trts_mix5,file='data/sex_trts_mix5.RData')

pb_upload("data/sex_trts_mix5.RData",
          name="sex_trts_mix5.RData",
          tag="data.v.1")


