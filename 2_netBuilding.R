## 2: Network (matrix) building from individual interactions.

## Loads the compiled dataset, generates simulated data randomizing sex,
## then breaks the observations into network matrices at the Site+Year level. 
## Lastly, it calculates node-level network metrics
## Modified (very slightly) from Ponisio sky islands data prep.


##this should be run on lauren's computer via terminal (git bash):
#ssh gsmith@osmia.dyn.ucr.edu
#ssh 192.168.1.5
#cd Documents
#git clone https://github.com/gosmith2/networks_by_sex.git
#cd Documents/networks_by_sex

## setwd("~/Dropbox/networks_by_sex")

pb_download("spec.all.RData",
            dest="data",
            tag="data.v.1")

load("data/spec.all.RData")

####-------------------------------------------####
#### Randomizing males and females for analyses
####-------------------------------------------####

## WARNING: the following code is computationally intesnive.
## It is currently written to be processed by 10 cores
## via the function mclapply. 

pb_download("spec.all.RData",
          dest="data",
          tag="data.v.1")

## cores should be 3 for bombus, 10 for osmia
cores <- 10

## randomize the sexes w/in species w/in sites. observed data is
## element 1 of the resulting list, randomized are elements 2-1000

rand_sexes <- ran.gen(spec.all, 999, cores)

save(rand_sexes,"data/rand_sexes.RData")



#build the networks at the sex level using the mixed sexes
nets_mix <- mclapply(rand_sexes, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)


## confirm that the networks are actually different
ifelse(any(nets_mix[[1]][[1]]!=nets_mix[[2]][[1]]),
       print("SUCCESS: the networks are different"),
       print("WARNING: the networks may not be different"))


#remove all networks with too few interactions to calculate metrics
nets_mix_clean <- mclapply(nets_mix, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)


#save the networks themselves and upload them
save(nets_mix_clean, file = 'data/nets_mix_clean.Rdata')
pb_upload("data/nets_mix_clean.RData",
          name="nets_mix_clean.RData",
          )

####-------------------------------------------####
#### Calculate node-level network parameters
####-------------------------------------------####

## calculate network stats at the node level, output into usable data frame
# reshuffling threshold is 5 in this case (i.e. there must be at least 5
# males and 5 females of a given species in a given network for it to be included)
sex_trts_mix5 <- mclapply(nets_mix_clean,
                          function(x) calcSpec(x, indiv = 5),
                          mc.cores = cores)


## confirm that the values are different
ifelse(any(sex.trts.mix5[[1]]$species.strength!=
             sex.trts.mix5[[2]]$species.strength),
       print("SUCCESS: the network statistics are different between randomizations"),
       print("WARNING: the network statistics may be not different between randomizations")
)


#save and upload
save(sex_trts_mix5,file='data/sex_trts_mix5.RData')

pb_upload("data/sex_trts_mix5.RData",
          name="sex_trts_mix5.RData",
          tag="data.v.1")



######-------------------------------

# - Alternative code that will be removed before pub

######-------------------------------

#rather than group observations at the Sp+site+yr level, grouping across sites w/in a year:
rand.sexes.SpYr.ls <- ran.gen(spec.all,3,cores,"SpYr")

nets.mix.SpYr <- mclapply(rand.sexes.SpYr.ls, function(y){
  breakNetMixSpYr(y, 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)

nets.mix.cleanSpYr <- mclapply(nets.mix.SpYr, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

sex.trts.mix2SpYr <- mclapply(nets.mix.cleanSpYr,
                              function(x) calcSpec(x, indiv = 2, lvl="SpYr"),
                              mc.cores = cores)


#Alternative reshuffling thresholds:
sex.trts.mix2 <- mclapply(nets.mix.clean,
                          function(x) calcSpec(x, indiv = 2),
                          mc.cores = cores)

sex.trts.mix3 <- mclapply(nets.mix.clean,
                          function(x) calcSpec(x, indiv = 3),
                          mc.cores = cores)

sex.trts.mix10 <- mclapply(nets.mix.clean,
                           function(x) calcSpec(x, indiv = 10),
                           mc.cores = cores)

sex.trts.mix20 <- mclapply(nets.mix.clean,
                           function(x) calcSpec(x, indiv = 20),
                           mc.cores = cores)

save(sex.trts.mix3,file='data/sex_trts_mix3.RData')
save(sex.trts.mix10,file='data/sex_trts_mix10.RData')
save(sex.trts.mix20,file='data/sex_trts_mix10.RData')



