#Network level analyses



####----------------------------####
#### Prep for networklevel analyses
####----------------------------####

pb_download("spec.all.RData",
            dest="data",
            tag="data.v.1")

load('data/spec.all.RData')

metric.net <- c('connectance',
                'number of compartments',
                'nestedness',
                'NODF',
                'H2',
                'robustness',
                'vulnerability',
                'niche overlap',
                'functional complementarity')

N = 999
cores = 10

##all networks, with males and females included separately
nets.obs.sex<-breakNetMix(spec.all,'Site','Year',"GenusSpeciesSex")

netStats.sex <- mclapply(nets.obs.sex,
                         calcNetworkMetrics,
                         N=N,
                         index=metric.net,
                         mc.cores=cores)

netStats.sexCI <- mclapply(nets.obs.sex,
                         calcNetworkMetricsCI,
                         N=N,
                         index=metric.net,
                         mc.cores=cores)


save(netStats.sex, file = "data/netStats_sex.RData")
save(netStats.sexCI, file = "data/netStats_sex_CI.RData")

pb_upload("data/netStats_sex.RData",
          name="netStats_sex.R=Data",
          tag="data.v.1")

pb_upload("data/netStats_sex_CI.RData",
          name="netStats_sex_CI.R=Data",
          tag="data.v.1")


##all networks, males and females lumped into species
nets.obs.sp <- breakNet(spec.all,'Site','Year')
netStats.sp <- mclapply(nets.obs.sp,
                        calcNetworkMetrics,
                        N=N,
                        index=metric.net,
                        mc.cores=cores)

netStats.spCI <- mclapply(nets.obs.sp,
                        calcNetworkMetricsCI,
                        N=N,
                        index=metric.net,
                        mc.cores=cores)


save(netStats.spCI, file = "data/netStats_sp_CI.RData")

pb_upload("data/netStats_sp.RData",
          name="netStats_sp.RData",
          tag="data.v.1")
pb_upload("data/netStats_sp_CI.RData",
          name="netStats_sp_CI.RData",
          tag="data.v.1")

##all networks, males dropped
spec.all.drop <- filter(spec.all, spec.all$Sex=="f")
nets.obs.f <- breakNet(spec.all.drop,'Site','Year')

netStats.f <- mclapply(nets.obs.f,
                       calcNetworkMetrics,
                       N=N,
                       index=metric.net,
                       mc.cores=cores)
netStats.fCI <- mclapply(nets.obs.f,
                       calcNetworkMetricsCI,
                       N=N,
                       index=metric.net,
                       mc.cores=cores)

save(netStats.f, file = "data/netStats_f.RData")
save(netStats.fCI, file = "data/netStats_f_CI.RData")

pb_upload("data/netStats_f.RData",
          name="netStats_f.RData",
          tag="data.v.1")
pb_upload("data/netStats_f_CI.RData",
          name="netStats_f_CI.RData",
          tag="data.v.1")



##Load and merge the 3 stats datasets
pb_download("netStats_f.RData",
            dest="data",
            tag="data.v.1")
pb_download("netStats_sex.RData",
            dest="data",
            tag="data.v.1")
pb_download("netStats_sp.RData",
            dest="data",
            tag="data.v.1")

load("data/netStats_sp.RData")
load("data/netStats_f.RData")
load("data/netStats_sex.RData")


names(netStats.f)

netStats.all <- do.call(rbind,
                        lapply(names(netStats.f),
                               function (x) {
                                 
                                 all.df <- data.frame(rbind(netStats.sex[x][[1]][1:39],
                                                            netStats.sp[x][[1]][1:39],
                                                            netStats.f[x][[1]][1:39]))
                                 names(all.df) <- names(netStats.sex["Zamora.2014"][[1]])
                                 all.df$SiteYr <- x
                                 all.df$trt <- c("sex","sp","fem")
                                 
                                 
                                 return(all.df)
                               })
)

ggplot(netStats.all, aes(x=pH2, color=trt)) +
  geom_density()

ggplot(netStats.all, aes(x=zniche.overlap.HL, color=trt)) +
  geom_density()

ggplot(netStats.all, aes(x=robustness.HL, color=trt)) +
  geom_density()

ggplot(netStats.all, aes(x=functional.complementarity.LL, color=trt)) +
  geom_density()

ggplot(netStats.all, aes(x=generality.HL, color=trt)) +
  geom_density()

ggplot(netStats.all, aes(x=zvulnerability.LL, color=trt)) +
  geom_density()


trtComp <- function (trt1,trt2) {
  one <- subset(netStats.all, netStats.all$trt == trt1)
  two <- subset(netStats.all, netStats.all$trt == trt2) 
  
  loop.ls <- c(1:39)
  combo <- lapply(loop.ls, function(x){
    comb <- as.vector(one[x]-two[x])
    return(comb)
  })
  
  clean <- do.call(cbind,combo)
  clean$SiteYr <- one$SiteYr
  clean$trt <- paste(trt1,trt2,sep = "-")
  return(clean)
}

sp_f.comp <- trtComp("sp","fem")
sp_sex.comp <- trtComp("sp","sex")





#______________________________________________________

#______________________________________________________




#######################################################

## Network level traits

#######################################################

pb_download("netlvlYH.RData",
            dest="data",
            tag="data.v.1")

load("data/netlvlYH.RData")

overallTest(netlvl,metric.net,zscore=F)


metric.net <- c('connectance',
                'number.of.compartments',
                'nestedness',
                'NODF',
                'robustness.HL',
                'robustness.LL',
                'vulnerability.LL',
                'H2',
                'niche.overlap.HL',
                'niche.overlap.LL',
                'functional.complementarity.HL',
                'functional.complementarity.LL')

calcNullPropNets <- function(data, metrics, zscore=TRUE) {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(dist.df$SiteYr), function(y) {
    #browser()
    site <- filter(dist.df, dist.df$SiteYr == y)
    obs <- filter(site, site$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      #browser()
      if (zscore == TRUE){
        metZ <- scale(site[,z],center = TRUE, scale = TRUE)
        #metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1])
        #gotta think about NA treatment here too. I kinda think its
        #getting rid of stuff again. Though this may be fixed when
        #I fix stuff above
      }else{
        metprop <- sum(site[,z] <= obs[,z]) / length(site$SiteYr)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SiteYr <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}

diffsNets.df <- calcNullPropNets(netlvl,metric.net,zscore=F) #quick

overallTest(diffsNets.df,metric.net,zscore=F)
#yes, def some differences here. not all of them really really sig, but for sure



genNullDistNets <- function(data, metrics, mean.by,zscore=TRUE) {
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    #browser()
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #extract an average value w/in each mean.by
  dist.build <- mclapply(unique(dist.df[,mean.by]), function(y) {
    net <- filter(dist.df, dist.df[,mean.by] == y)
    obs <- filter(net, net$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      browser()
      if (zscore == TRUE){
        mean.Z <- mean(scale(net[,z],center=T,scale=T),na.rm=T)
      }else{
        mean.value <- mean(net[,z])
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$mean.by <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,dist.build)
  return(sig.dist)
}

nullDistNets.df <- genNullDistNets(netlvl,metric.net,"sim",zscore=F)

save(nullDistNets.df,file="data/nullDistNetsYH.RData")

pb_upload("data/nullDistNetsYH.RData",
          name="nullDistNetsYH.RData",
          tag="data.v.1")
pb_download("nullDistNetsYH.RData",
            dest="data",
            tag="data.v.1")

meanObsDiffNet<-lapply(metric.net, function(x){
  mean(netlvl[[1]][,x],na.rm=T)
})

plot(density(nullDistNets.df$robustness.HL,na.rm = T))
abline(v=meanObsDiffNet[5])








####################################

##Below was exploratory, not run now

####################################


#bargraph.CI(response=sex_trts.df$d,x.factor=sex_trts.df$sex)

#Exploratory:
##females higher degree (1.8 vs 1.45ish). same pattern diff #s for normalized
##**females higher str (0.4 vs. 0.25)
##males lower push pull (~-.7 vs -.6)
##males slightly higher nestedrank (0.55 vs 0.47ish)
##PDI about the same, v little var in either
##males slightly higher for resource.range (.92 vs .89? ish)
##sp. spec index about as above
##females higher PSI (0.19 vs .15)
##males v slightly higher node spec index (0.15ish vs 0.155?)
##**females much more between than males (0.055 vs 0.025). weighted even more diff
##females slightly more close than males (0.05 vs 0.047?). larger diff for weighted
##males higher fisher A (114 mil vs 96 mil)
##**females higher partner diversity (~0.3 vs 0.2)
##females higher eff. partners (1.58 vs 1.38)
##females slightly higher prop generality (0.4 vs 0.35)
##females slightly higher prop similarity (0.42 vs 0.38)
## d essentially the same





###--------------
## Generate null distributions and plotting

#this one averages everything by simulation
nullDist.df <- genNullDist(sexDiffs2.df,metric.ls,
                           zscore=T)
save(nullDist.df,file='data/nullDist.RData')

#averages by sp+site+year (how overallTest was run)
nullDistDiff2.df <- genNullDist(sexDiffs2.df,metric.ls,
                                zscore=F)

save(nullDistDiff.df,file='data/nullDistDiff.RData')


pb_upload("data/nullDist.RData",
          name="nullDist.RData",
          tag="data.v.1")

pb_download("nullDist.RData",
            dest="data",
            tag="data.v.1")

pb_upload("data/nullDistDiff.RData",
          name="nullDistDiff.RData",
          tag="data.v.1")

pb_download("nullDistDiff.RData",
            dest="data",
            tag="data.v.1")