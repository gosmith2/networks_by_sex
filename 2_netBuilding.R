## Network (matrix) building from individual interactions.
## Modified (very slightly) from Ponisio sky islands data prep.


##this should be run on lauren's computer via terminal (git bash):
#ssh gsmith@osmia.dyn.ucr.edu
#cd Documents
#git clone https://github.com/gosmith2/networks_by_sex.git
#cd Documents/networks_by_sex

## setwd("~/Dropbox/networks_by_sex")

pb_download("spec.all.RData",
            dest="data",
            tag="data.v.1")

load("data/spec.all.RData")

####-------------------------------------------####
#### Randomizing males and females for specieslevel analyses
####-------------------------------------------####

## WARNING: the following code is computationally intesnive.
## It is currently written to be processed by 10 processor
## cores via the function mclapply. 

pb_download("spec.all.RData",
          dest="data",
          tag="data.v.1")

## cores should be 3 for bombus, 10 for osmia
cores <- 10

## randomize the sexes w/in species w/in sites. observed network is
## element 1 of the resulting list, randomized are elements 2-1000

rand.sexes.ls <- ran.gen(spec.all, 999, cores)


#build the networks at the sex level using the mixed sexes
nets.mix <- mclapply(rand.sexes.ls, function(y){
  breakNetMix(y, 'Site', 'Year', 'GenusSpeciesMix')
}, mc.cores = cores)

## confirm that the networks are actually different
ifelse(nets.mix[[1]]$Zamora.2014[[8]]
       ==
       nets.mix[[2]]$Zamora.2014[[8]],
       print("WARNING: the networks may be not different between randomizations"),
       print("SUCCESS: the networks are different between randomizations")
)

#remove all networks with too few interactions to calculate metrics
nets.mix.clean <- mclapply(nets.mix, function(x){
  x[sapply(x, function(y) all(dim(y) > 1))]
}, mc.cores=cores)

#save the networks themselves
save(nets.mix.clean, file = 'data/mix_netsYHS.RData')

####-------------------------------------------####
#### Calculate node-level network parameters
####-------------------------------------------####

#calculate network stats at the node level, output into usable data frame
sex.trts.mix2 <- mclapply(nets.mix.clean,
                          function(x) calcSpec(x, indiv = 2),
                          mc.cores = cores)

## confirm that the values are different
ifelse((sex.trts.mix2[[1]] %>%
	filter(Site == "Zamora", Year == 2014) %>%
  select(species.strength) %>%
  head())
	==
	  (sex.trts.mix2[[2]] %>%
	     filter(Site == "Zamora", Year == 2014) %>%
	     select(species.strength) %>%
	     head()),
	print("WARNING: the network statistics may be not different between randomizations"),
	print("SUCCESS: the network statistics are different between randomizations")
)
	

save(sex.trts.mix2,file='data/sex_trts_mix2.RData')


pb_upload("data/sex_trts_mixYHS2.RData",
          name="sex_trts_mixYHS2.RData",
          tag="data.v.1")
pb_upload("data/mix_netsYHS.RData",
			name="mix_netsYHS.RData",
            tag="data.v.1")


