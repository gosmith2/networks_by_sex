## 2: Network (matrix) building from individual interactions.

## Loads the compiled dataset, generates simulated data randomizing sex,
## then breaks the observations into network matrices at the Site+Year level. 
## Lastly, it calculates node-level network metrics
## Modified (very slightly) from Ponisio sky islands data prep.

## all working directories are relative to the gitrepo 
library(parallel)
source("prepNets.R")

load("data/spec_all.RData")

####-------------------------------------------####
#### Randomizing males and females for analyses
####-------------------------------------------####

## WARNING: the following code is computationally intesnive.
## It is currently written to be processed by 10 cores
## via the function mclapply. 


## the number of cores you wish to run the simulation in parallel across 
cores <- 10

## From the observation data, generate a validation dataframe ("S" for "same") where males and
## females share the same plant visitation vector. "O" (for "original") is also run through
## this process, but visitation is left unchanged
## O and S will be run through the proceeding steps in parallel. 
spec_allS <- plant.shuffler(spec_all,trt='same')
spec_allO <- plant.shuffler(spec_all,trt='orig')


## randomize the sexes w/in species w/in sites. observed data is
## element 1 of the resulting list, randomized are elements 2-1000
rand_sexes2kS1 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kS2 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kS3 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kS4 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kS5 <- ran.gen(spec_allS, 2000, cores)

rand_sexes2kO1 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kO2 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kO3 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kO4 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kO5 <- ran.gen(spec_allO, 2000, cores)

##Note: 2000 was the highest number of iterations that successfully compiled as a single object, 
##presumably due to memory constraints. Multiple objects were therefore created to get to 
##higher iteration numbers


#save the objects to pause and load later
save(rand_sexes2kO1,file='data/rand_sexes2kO1.RData')
save(rand_sexes2kO2,file='data/rand_sexes2kO2.RData')
save(rand_sexes2kO3,file='data/rand_sexes2kO3.RData')
save(rand_sexes2kO4,file='data/rand_sexes2kO4.RData')
save(rand_sexes2kO5,file='data/rand_sexes2kO5.RData')

save(rand_sexes2kS1,file='data/rand_sexes2kS1.RData')
save(rand_sexes2kS2,file='data/rand_sexes2kS2.RData')
save(rand_sexes2kS3,file='data/rand_sexes2kS3.RData')
save(rand_sexes2kS4,file='data/rand_sexes2kS4.RData')
save(rand_sexes2kS5,file='data/rand_sexes2kS5.RData')



load('data/rand_sexes2kS1.RData')
load('data/rand_sexes2kS2.RData')
load('data/rand_sexes2kS3.RData')
load('data/rand_sexes2kS4.RData')
load('data/rand_sexes2kS5.RData')



#build the networks at the sex level using the mixed sexes
nets_mix_S1 <- mclapply(rand_sexes2kS1, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S2 <- mclapply(rand_sexes2kS2, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S3 <- mclapply(rand_sexes2kS3, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S4 <- mclapply(rand_sexes2kS4, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S5 <- mclapply(rand_sexes2kS5, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)

nets_mix_O1 <- mclapply(rand_sexes2kO1, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O2 <- mclapply(rand_sexes2kO2, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O3 <- mclapply(rand_sexes2kO3, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O4 <- mclapply(rand_sexes2kO4, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O5 <- mclapply(rand_sexes2kO5, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)


## confirm that the networks are actually different
ifelse(any(nets_mix_S1[[1]][[1]]!=nets_mix_S1[[2]][[1]]),
       print("SUCCESS: the networks are different"),
       print("WARNING: the networks are not different"))


##if multiple objects were required to get to desired simulation number, this combines them. 
##The network objects are smaller and don't cause as many memory issues when combined. 
##NOTE: element 1 in all of the lists is the seed data (i.e., either Original or Same), so was
## only included in the first list object. 
nets_mix_S10k<- c(nets_mix_S1,
                      nets_mix_S2[2:2001],
                      nets_mix_S3[2:2001],
                      nets_mix_S4[2:2001],
                      nets_mix_S5[2:2001])

nets_mix_O10k<- c(nets_mix_O1,
                  nets_mix_O2[2:2001],
                  nets_mix_O3[2:2001],
                  nets_mix_O4[2:2001],
                  nets_mix_O5[2:2001])


#remove all networks with too few interactions to calculate metrics
nets_mix_clean10kO <- mclapply(nets_mix_O10k, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

nets_mix_clean10kS <- mclapply(nets_mix_S10k, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

#save the networks
save(nets_mix_clean10kO, file = 'data/nets_mix_clean10kO.RData')

save(nets_mix_clean10kS, file = 'data/nets_mix_clean10kS.RData')


#save smaller lists of the networks for use later / github
  #Note: even at 500 reps these files are too large for github
nets_mix_clean_SmallO <- nets_mix_clean10kO[1:400]
save(nets_mix_clean_SmallO,file='data/nets_mix_clean_SmallO.RData')

nets_mix_clean_SmallS <- nets_mix_clean10kS[1:400]
save(nets_mix_clean_SmallS,file='data/nets_mix_clean_SmallS.RData')

