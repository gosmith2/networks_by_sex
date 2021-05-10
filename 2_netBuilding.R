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



####BELOW IS THE STEP THAT NOW DOES THE BOOTSTRAPPING
  #key for bootstraps already done (with rand_sexes <- "x")
  # x2k = no bootsrapping, original
  # x1kboot[number] = all sites bootstrapped by the same number. 
    #NOTE: x1kboot10 was the first attempt; some older objects may still be named just "x1kboot" rather than "x1kboot10"
    #in current code, this would be replicated by havin the bootnum vector be 3 of the same number (e.g., c(5,5,5))
  # x2kbootST = stood for Standardized Time, but was actually bootnum=c(2,3,1) for yos, hr, and si
  # x2kbootT = stands for Time, and is bootnum = c(4,6,1)
  # these strings after x should carry through to all downstream objects. 


rand_sexes1kbootT <- ran.gen(spec_all, 1000, cores,boot=TRUE,bootnum=c(4,6,1))
  #Note: 2000 was the highest number of iterations that successfully compiled, 
  #presumably due to memory constraints. 



save(rand_sexes1kbootST,file="data/rand_sexes1kbootST.RData")

#build the networks at the sex level using the mixed sexes
nets_mix_bootT4 <- mclapply(rand_sexes2kbootT4, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)


## confirm that the networks are actually different
ifelse(any(nets_mix_bootST[[1]][[1]]!=nets_mix_bootST[[2]][[1]]),
       print("SUCCESS: the networks are different"),
       print("WARNING: the networks are not different"))


##if multiple objects were required to get to 1k or 2k simulations, this combined them
#nets_mix_2kbootT <- c(nets_mix_bootT1,
#                      nets_mix_bootT2[2:501],
#                      nets_mix_bootT3[2:501],
#                      nets_mix_bootT4[2:501])
                          #  nets_mix_boot203,
                          #  nets_mix_boot204,
                          #  nets_mix_boot205,
                          #  nets_mix_boot206,
                          #  nets_mix_boot207)

#remove all networks with too few interactions to calculate metrics
nets_mix_clean2kbootT <- mclapply(nets_mix_2kbootT, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

#save the networks themselves and upload them
save(nets_mix_clean2kbootT, file = 'data/nets_mix_clean2kbootT.RData')
pb_upload("data/nets_mix_clean2kbootT.RData",
          name="nets_mix_clean2kbootT.RData",
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
sex_trts_mix5bootT <- mclapply(nets_mix_clean2kbootT,
                          function(x) calcSpec(x, indiv = 5),
                          mc.cores = cores)

## confirm that the values are different
ifelse(any(sex_trts_mix5bootST[[1]]$weighted.closeness!=
             sex_trts_mix5bootST[[2]]$weighted.closeness),
       print("SUCCESS: the network statistics are different between randomizations"),
       print("WARNING: the network statistics are not different between randomizations")
)

#save and upload
save(sex_trts_mix5bootT,file='data/sex_trts_mix5bootT.RData')

pb_upload("data/sex_trts_mix5bootST.RData",
          name="sex_trts_mix5bootST.RData",
          tag="data.v.1")


