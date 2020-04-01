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


##### NOTE!!!!! Before running these, you should probably subset down to the networks where at least 1 species 
  # displayed sex differences!! e.g.:

#spec.net <- spec.all[spec.all$site %in% unique(list)]

#AND REMEMBER TO CHANGE SPEC.ALL to SPEC.NET DOWNSTREAM!!


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
          name="netStats_sex_CI.RData",
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


#######################

#load CI versions 

#######################

pb_download("netStats_f_CI.RData",
            dest="data",
            tag="data.v.1")
pb_download("netStats_sex_CI.RData",
            dest="data",
            tag="data.v.1")
pb_download("netStats_sp_CI.RData",
            dest="data",
            tag="data.v.1")
load("data/netStats_sp_CI.RData")
load("data/netStats_f_CI.RData")
load("data/netStats_sex_CI.RData")



netStatsCI.all <- do.call(rbind,
                        lapply(names(netStats.fCI)[-c(116,179)],
                               function (x) {
                                 
                                # browser()
                                 
                                 if(any(!is.na(netStatsCI.sex[[x]]))) {
                                   fullSex.df <- as.data.frame(t(netStatsCI.sex[[x]]))
                                   fullSex.df$site <- x
                                   fullSex.df$trt <- "sex" 
                                 } else {
                                   fulls.df <- as.data.frame(t(netStatsCI.sex[[x]]))
                                   rownames(fullSex.df) <- rownames(t(netStatsCI.sex[[1]]))
                                   colnames(fullSex.df) <- colnames(t(netStatsCI.sex[[1]]))
                                   fullSex.df$site <- x
                                   fullSex.df$trt <- "sex"
                                 }
                                 
                                 if(any(!is.na(netStats.fCI[[x]]))) {
                                   fullf.df <- as.data.frame(t(netStats.fCI[[x]]))
                                   fullf.df$site <- x
                                   fullf.df$trt <- "f" 
                                 } else {
                                   fullf.df <- as.data.frame(t(netStats.fCI[[x]]))
                                   rownames(fullf.df) <- rownames(t(netStatsCI.sex[[1]]))
                                   colnames(fullf.df) <- colnames(t(netStatsCI.sex[[1]]))
                                   fullf.df$site <- x
                                   fullf.df$trt <- "f"
                                 }
                                 
                                 if(any(!is.na(netStats.spCI[[x]]))) {
                                   fullsp.df <- as.data.frame(t(netStats.spCI[[x]]))
                                   fullsp.df$site <- x
                                   fullsp.df$trt <- "sp" 
                                 } else {
                                   fullsp.df <- as.data.frame(t(netStats.spCI[[x]]))
                                   rownames(fullsp.df) <- rownames(t(netStatsCI.sex[[1]]))
                                   colnames(fullsp.df) <- colnames(t(netStatsCI.sex[[1]]))
                                   fullsp.df$site <- x
                                   fullsp.df$trt <- "sp"
                                 }
                                 
                                 
                                 all.df <- data.frame(rbind(fullSex.df,
                                                            fullf.df,
                                                            fullsp.df))
                                 all.df$value <- rownames(all.df)
                                 rownames(all.df) <- NULL
                                 all.df$value <- rep(c("95L","95H","90L","90H","obs"))
                                 return(all.df)
                               })
)


netsize.df <- as.data.frame(unique(netStatsCI.all$site))
netsize.df$plants <- unlist(lapply(unique(netStatsCI.all$site), function(x){
  dimensions <- dim(nets.obs.sp[[x]])
  return(dimensions[1])
  }))
netsize.df$total <- unlist(lapply(unique(netStatsCI.all$site), function(x){
  dimensions <- dim(nets.obs.sp[[x]])
  return(sum(dimensions))
}))
names(netsize.df) <- c("site","plants","total")


### add netsize columns to netStatsCI.all using matching
netStats <- left_join(netStatsCI.all, netsize.df, by="site")
names(netStats) <- c(names(netStats)[1:13],"network",names(netStats)[15:18])
netStats <- mutate(netStats, site = sub("\\.[0-9]+$", "", network))

#subset down to the sites where at least one species showed 
sexnetstats <- netStats[netStats$network %in% names(nets.mix.clean[[1]]),]

### plot values (e.g., robustness) vs network size, look for interactions

ggplot(sexnetstats)+
#  geom_line(aes(x=robustness.HL, y=plants, color=trt),data= .%>% filter(value=="obs")) +
  geom_point(aes(x=robustness.HL,y=site,color=trt),data= .%>% filter(value=="obs"))
  #hard to see whats up


### "test" mains as proportions of sites outside of CIs


### also test by just regressing the obs robustness (etc) by trt to look at whether there are overall trends
  #could occur independently of whether too many are "significantly" diff from 

trt_rob.lm <- lm(robustness.HL~trt*plants,data=sexnetstats, subset = sexnetstats$value=="obs")
summary(trt_rob.lm)
  #sig; sex and sp seem to be larger than f
  #no real strong effect of network size. very very slight negative relationship (which is sig though)
  
  #posthoc for sex vs sp
trt_rob_sexsp.lm <- lm(robustness.HL~trt*plants,data=sexnetstats, subset = sexnetstats$value=="obs" & sexnetstats$trt!="f")
summary(trt_rob_sexsp.lm)
  #sp is marginally sig (0.09) smaller than sex. very little diff though

#lmm w/ site
trt_rob_site.lmm <- lme(robustness.HL~trt,data=sexnetstats, subset = sexnetstats$value=="obs", random=~1|site)
summary(trt_rob_site.lmm)
  #like... exactly the same as above

#lmm w/ size?  
trt_rob_size.lmm <- lme(robustness.HL~trt,data=sexnetstats, subset = sexnetstats$value=="obs", random=~1|plants)
summary(trt_rob_size.lmm)
  #also exactly the same as above. huh....?



### calculate the proportions diff for each metric: 

CItest <- function(data, CI = "90", comp){
  nets <- lapply(unique(data$network),function(x){
    
    #browser()
    #subset to the network
    net <- filter(sexnetstats, network == x)
    
    #pull out the CI values we'll compare
    spCompH <- filter(net,trt=="sp",value==paste0(CI,"H"))
    spCompL <- filter(net,trt=="sp",value==paste0(CI,"L"))
    compH <- filter(net, trt==comp,value==paste0(CI,"H"))
    compL <- filter(net, trt==comp,value==paste0(CI,"L"))
    
    #make comparisons for each statistic
    stats <- lapply(names(sexnetstats)[1:13],function(y){
      higher <- spCompH[,y] - compL[,y] #only way this is negative is if the comparison is significantly larger
      lower <- compH[,y] - spCompL[,y] #only way this is negative is if the comparison is significantly smaller
      
      stat.df <- as.data.frame(higher)
      stat.df$lower <- lower
      names(stat.df) <- c(paste(y,"higher",sep="_"),paste(y,"lower",sep="_"))
      return(stat.df)
    })
    
    stats <- do.call(cbind,stats)
    stats$network <- x
    return(stats)
  })
  
  nets <- do.call(rbind, nets)
  return(nets)

}

test95 <- CItest(sexnetstats, CI="95",comp="sex")
sig.vec95 <- unlist(lapply(names(test95)[1:26],function(x){
  sum(test95[,x]<0)/length(test95[,x])
}))
names(sig.vec95) <- names(test95)[1:26]
sig.vec95

#Woa, some interesting things in here. value = proportion of networks for which that value (e.g., connectance) was
#either higher for the sex lvl networks (relative to sp lvl) or lower. 

  #connectance was much more frequently lower for sex networks (2.7% higher vs 97% lower)
    #Makes sense: you're adding new pol nodes, but keeping number of interactions the same

  #more compartments in sex networks (65% higher, 0.7% lower). 
    #also makes sense for same reason as above

  #no strong trend for nestedness (47% higher, 28% lower)

  #NODF was much more frequently lower for sex networks (3.5% higher, 95% lower)
    #??

  #H2 more frequently higher (81% vs 7%)

  #pollinator niche overlap was lower (3.5% vs 70% lower)
    #again makes sense: you're dividing interactions into smaller sub-nodes, -> by numbers fewer "generalists" w/ overlap

  #same story for plant niche overlap (0.4% vs 90% lower)

  #Pollinator robustness was more often higher in sex networks (70% of networks showed higher vs 11$ lower)
    #not sure exactly why

  #plant robustness was much LOWER in sex networks (4.7% higher vs 93% of networks where it was lower)
    #gotta think about why

  #fxn comp for polls: not too strong but seems lower: 20% higher vs 63$ lower

  #plant fxnal comp was lower: 12% higher, 79% lower

  #generality was lower for sex networks: 1.6% higher vs 96% lower
    #again, makes sense. splitting -> more "specialists"

  #plant vulnerability was definitely higher (97% vs 0.0000000%)


### test interactions by regressing those proportions against size

test95f <- CItest(sexnetstats, CI="95",comp="f")
sig.vec95f <- unlist(lapply(names(test95f)[1:26],function(x){
  sum(test95f[,x]<0)/length(test95f[,x])
}))
names(sig.vec95f) <- names(test95f)[1:26]
sig.vec95f




myplot <- ggplot(data=netStatsCI.all)+geom_density(aes(x=robustness.HL, color=trt))
myplot %+% subset(netStatsCI.all, value %in% c("obs"))
myplot %+% subset(netStatsCI.all, value %in% c("95H"))

ggplot(data=netStatsCI.all)+
  geom_density(aes(x=robustness.HL, color=trt),data= . %>% filter(value=="obs")) +
  geom_density(aes(x=robustness.HL, color=trt,linetype="dashed"),data= . %>% filter(value=="95H")) 


geom_line(aes(Value1, Value2, group=ID, colour=ID),
          ,subset = .(ID %in% c("P1" , "P3")))

myplot<-ggplot(df)+geom_line(aes(Value1, Value2, group=ID, colour=ID))
myplot %+% subset(df, ID %in% c("P1","P3"))
myplot %+% subset(df, ID %in% c("P2"))

ggplot(netStats.all, aes(x=zvulnerability.LL, color=trt)) +
  geom_density()


ggplot(data=netStatsCI.all)+
  geom_density(aes(x=robustness.HL, color=trt),data= . %>% filter(value=="obs")) +
  geom_density(aes(x=robustness.HL, color=trt,linetype="dashed"),data= . %>% filter(value=="95H")) 


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