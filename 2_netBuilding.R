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

## From the observation data, generate a validation dataframe ("S" for "same") where males and
## females share the same plant visitation vector. "O" (for "original") is also run through
## this process, but visitation is left unchanged
## O and S will be run through the proceeding steps in parallel. 
spec_allS <- plant.shuffler(spec_all,trt='same')
spec_allO <- plant.shuffler(spec_all,trt='orig')


## randomize the sexes w/in species w/in sites. observed data is
## element 1 of the resulting list, randomized are elements 2-1000
rand_sexes2kNEWS1 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kNEWS2 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kNEWS3 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kNEWS4 <- ran.gen(spec_allS, 2000, cores)
rand_sexes2kNEWS5 <- ran.gen(spec_allS, 2000, cores)

rand_sexes2kNEWO1 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kNEWO2 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kNEWO3 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kNEWO4 <- ran.gen(spec_allO, 2000, cores)
rand_sexes2kNEWO5 <- ran.gen(spec_allO, 2000, cores)

##Note: 2000 was the highest number of iterations that successfully compiled as a single object, 
##presumably due to memory constraints. Multiple objects were therefore created to get to 
##higher iteration amounts

#rand_sexes1kNEWO <- ran.gen(spec_allO, 1000, cores)#,boot=TRUE,bootnum=c(4,6,1))


save(rand_sexes2kNEWO1,file='data/rand_sexes2kNEWO1.RData')
save(rand_sexes2kNEWO2,file='data/rand_sexes2kNEWO2.RData')
save(rand_sexes2kNEWO3,file='data/rand_sexes2kNEWO3.RData')
save(rand_sexes2kNEWO4,file='data/rand_sexes2kNEWO4.RData')
save(rand_sexes2kNEWO5,file='data/rand_sexes2kNEWO5.RData')

save(rand_sexes2kNEWS1,file='data/rand_sexes2kNEWS1.RData')
save(rand_sexes2kNEWS2,file='data/rand_sexes2kNEWS2.RData')
save(rand_sexes2kNEWS3,file='data/rand_sexes2kNEWS3.RData')
save(rand_sexes2kNEWS4,file='data/rand_sexes2kNEWS4.RData')
save(rand_sexes2kNEWS5,file='data/rand_sexes2kNEWS5.RData')



load('data/rand_sexes2kNEWS1.RData')
load('data/rand_sexes2kNEWS2.RData')
load('data/rand_sexes2kNEWS3.RData')
load('data/rand_sexes2kNEWS4.RData')
load('data/rand_sexes2kNEWS5.RData')



#build the networks at the sex level using the mixed sexes
nets_mix_S1 <- mclapply(rand_sexes2kNEWS1, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S2 <- mclapply(rand_sexes2kNEWS2, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S3 <- mclapply(rand_sexes2kNEWS3, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S4 <- mclapply(rand_sexes2kNEWS4, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_S5 <- mclapply(rand_sexes2kNEWS5, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)

nets_mix_O1 <- mclapply(rand_sexes2kNEWO1, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O2 <- mclapply(rand_sexes2kNEWO2, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O3 <- mclapply(rand_sexes2kNEWO3, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O4 <- mclapply(rand_sexes2kNEWO4, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)
nets_mix_O5 <- mclapply(rand_sexes2kNEWO5, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)


## confirm that the networks are actually different
ifelse(any(nets_mix_NEWS1[[1]][[1]]!=nets_mix_NEWS1[[2]][[1]]),
       print("SUCCESS: the networks are different"),
       print("WARNING: the networks are not different"))


##if multiple objects were required to get to desired simulation number, this combines them. 
##The network objects are smaller and don't cause as many memory issues
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

#save the networks themselves and upload them
save(nets_mix_clean10kO, file = 'data/nets_mix_clean10kO.RData')
pb_upload("data/nets_mix_clean10kO.RData",
          name="nets_mix_clean10kO.RData",
          )
save(nets_mix_clean10kS, file = 'data/nets_mix_clean10kS.RData')
pb_upload("data/nets_mix_clean10kS.RData",
          name="nets_mix_clean10kS.RData",
)
pb_download('nets_mix_clean10kS.RData',
            dest='data',
            tag="data.v.1")




