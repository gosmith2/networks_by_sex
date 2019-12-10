## the purpose of this function is to break up data with many
## sites/years and prepare it for network analysis.

dropNet <- function(z){
  z[!sapply(z, FUN=function(q){
    any(dim(q) < 3)
  })]
}



breakNet <- function(spec.dat, site, year){
  ## puts data together in a list and removes empty matrices
  agg.spec <- aggregate(list(abund=spec.dat$GenusSpecies),
                        list(GenusSpecies=spec.dat$GenusSpecies,
                             Site=spec.dat[,site],
                             Year=spec.dat[,year],
                             PlantGenusSpecies=
                               spec.dat$PlantGenusSpecies),
                        length)
  sites <- split(agg.spec, agg.spec[,site])
  networks <- lapply(sites, function(x){
    split(x, f=x[,"Year"])
  })
  ## formats data matrices appropriate for network analysis
  comms <- lapply(unlist(networks, recursive=FALSE), function(y){
    samp2site.spp(site=y[,"PlantGenusSpecies"],
                  spp=y[,"GenusSpecies"],
                  abund=y[,"abund"])
  })
  return(comms)
}

breakNetMix <- function(spec.dat, site, year, mix){
  ## puts data together in a list and removes empty matrices
  agg.spec <- aggregate(list(abund=spec.dat[,mix]),
                        list(GenusSpeciesMix=spec.dat[,mix],
                             Site=spec.dat[,site],
                             Year=spec.dat[,year],
                             PlantGenusSpecies=
                               spec.dat$PlantGenusSpecies),
                        length)
  sites <- split(agg.spec, agg.spec[,site])
  networks <- lapply(sites, function(x){
    split(x, f=x[,"Year"])
  })

  ## formats data matrices appropriate for network analysis
  comms <- lapply(unlist(networks, recursive=FALSE), function(y){
    samp2site.spp(site=y[,"PlantGenusSpecies"],
                  spp=y[,"GenusSpeciesMix"],
                  abund=y[,"abund"])
  })
  return(comms)
}


getSpecies <- function(networks, FUN){
  species.site <- lapply(networks, FUN)
  site.plant <- rep(names(species.site), lapply(species.site, length))
  species <- data.frame(species=do.call(c, species.site),
                        siteStatus=site.plant,
                        site= sapply(strsplit(site.plant, "_"),
                                     function(x) x[1]),
                        status= sapply(strsplit(site.plant, "_"),
                                       function(x) x[2]))
  return(species)
}



## calculates the degree of species in a network
getDegree <- function(x, MARGIN){
  apply(x, MARGIN, function(y){
    length(y[y != 0])/length(y)
  })
}

## calculate various stats
calcStats <- function(x){
  means=mean(x)
  medians=median(x)
  mins <- min(x)
  maxs <- max(x)
  sds <- sd(x)
  return(c(mean=means,
           median=medians,
           min=mins,
           max=maxs,
           sd=sds))
}



## number of species that interact
getCon <- function(x, INDEX){
  apply(x, INDEX, function(y) sum(y > 0))
}

## calculates the species roles from a network and returns a dataframe
## with site status and ypr
## takes networks and specimen data

calcSpec <- function(nets, spec, indiv = 1){
  ## applies specieslevel from bipartite to networks
  species.lev <- lapply(nets, function(x){
    
    #calculate values
    sl <- specieslevel(x)
    
    #keep only species with observations >= indiv
    x <- as.data.frame(x)
    sums <- lapply(x, sum)
    cullNames <- names(sums[sums >= indiv])

    sl$`higher level` <- sl$`higher level`[cullNames,]
    
    return(sl)
  })
  
  #sum up under species (rowsums or such)
  
  ## extract the values and make a dataframe
  specs  <-  mapply(function(a, b)
    getSpec(species.lev = a,
            names.net = b,
            seps="[.]"),
    a = species.lev,
    b = names(nets),
    SIMPLIFY = FALSE)
  
  specs <- do.call(rbind, specs)
  specs$sex <- substr(specs$GenusSpecies,
                      nchar(specs$GenusSpecies),
                      nchar(specs$GenusSpecies)
                      )
  specs$sex <- ifelse(specs$speciesType=="plant"|specs$sex=="e",NA,specs$sex)
  specs$sex <- as.factor(specs$sex)
  rownames(specs) <- NULL
  return(specs)
}


## extreact specialization scores from specieslevel function and
## return data frame
getSpec <- function(species.lev, names.net, seps="_"){
  n.pp <- sapply(species.lev, nrow)
  pp <- c(unlist(sapply(species.lev, rownames)))
  names(pp) <- NULL
  all.pp <- do.call(rbind, species.lev)
  rownames(all.pp) <- NULL
  try(all.pp$GenusSpecies <- pp)
  all.pp$speciesType <- c(rep("pollinator", n.pp[1]),
                          rep("plant", n.pp[2]))
  all.pp$Site <- strsplit(names.net, seps)[[1]][1]
  all.pp$Year <- strsplit(names.net, seps)[[1]][2]
  return(all.pp)
}


#build a large heirarcical list of networks where sex is randomized


removeNets <- function(spec.data) {
  #narrows down to each site+year combo, then only keeps the networks
  #where at least 1 species has both sexes present
  col.ls <- lapply(unique(spec.data$SiteYr),function(x){
    if(length(unique(spec.data$GenusSpeciesSex[spec.data$SiteYr==x]))
       >
       length(unique(spec.data$GenusSpecies[spec.data$SiteYr==x]))){
      filter(spec.data,SiteYr==x)
  }})
  col.ls <- col.ls[!sapply(col.ls, is.null)]
  return(do.call(rbind, col.ls))
}

#this one without the network remove code that I put into removeNets
ran.sex <- function(spec.data){
  col.ls <- lapply(unique(spec.data$SiteYr),function(x){
    net <- filter(spec.data,SiteYr==x)
    sp.mix.col <-lapply(unique(net$GenusSpecies), function(y){
      sp<-filter(net,net$GenusSpecies==y)
      if(length(unique(sp$Sex)) != 1){
        sp$MixSex<-sample(sp$Sex,replace=FALSE)
      }else{
        sp$MixSex <- sp$Sex
      }
      return(sp)
    })
  })
  return(do.call(rbind, unlist(col.ls, recursive=FALSE)))
}


ran.gen<-function(spec.data,iterations,cores){
  #setup: add column to actual observation df
  spec.data$MixSex<-spec.data$Sex
  spec.data$GenusSpeciesMix<-paste(spec.data$GenusSpecies,
                                   spec.data$MixSex,
                                   sep="_")
  
  #Keep only networks where at least one sp has both sexes
  spec.keeps <- removeNets(spec.data)
  
  #innitiate list
  destList<-list() 
  destList[[1]]<-spec.keeps
  
  #scramble sex column within each species within each network
  it.vec<-1:iterations
  randoms<-mclapply(it.vec, function(z){
    spec.keeps<-ran.sex(spec.keeps)
    spec.keeps$GenusSpeciesMix<-paste(spec.keeps$GenusSpecies,
                                     spec.keeps$MixSex,
                                     sep="_")
    return(spec.keeps)
  }, mc.cores = cores)
  destList<-c(destList,randoms)
  return(destList)
}


calcMetric <- function(dat.web, index) {
  ## calculates modularity
  dat.web <- as.matrix(empty(dat.web))
  mets <-  networklevel(dat.web, index=index)
  return(mets)
}


## function to simulate 1 null, and calculate statistics on it
calcNullStat <- function(dat.web,
                         null.fun,...) {
  sim.web <- null.fun(dat.web)
  return(calcMetric(sim.web,...))
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
calcNetworkMetrics <- function (dat.web, N,
                                index= c("mean number of links",
                                         "niche overlap",
                                         "fc")) {
  ## calculate pvalues
  pvals <- function(stats, nnull){
    rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
  }
  ## calculate zvalues two different ways
  zvals <-function(stats){
    z.sd <- (stats[,1] -
               apply(stats, 1, mean, na.rm = TRUE))/
      apply(stats, 1, sd, na.rm = TRUE)
    z.sd[is.infinite(z.sd)] <- NA
    return(z.sd)
  }
  ## check that matrix is proper format (no empty row/col and no NAs)
  if(all(is.na(dat.web) == FALSE)) {
    ## drop empty rows and columns
    dat.web <- as.matrix(empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(is.matrix(dat.web)){
      if(all(dim(dat.web) >= 2)) {
        ## calculate null metrics
        null.stat <- replicate(N,
                               calcNullStat(dat.web,
                                            null.fun= vaznull.fast,
                                            index=index),
                               simplify=TRUE)
        ## calculate metrics from data
        true.stat <- calcMetric(dat.web,
                                index=index)
        out.mets <- cbind(true.stat, null.stat)
        ## compute z scores
        zvalues <- zvals(out.mets)
        names(zvalues) <- paste("z", names(true.stat), sep="")
        ## compute p-values
        pvalues <- pvals(out.mets, N)
        names(pvalues) <- paste("p", names(true.stat), sep="")
        out <- c(true.stat, zvalues, pvalues)
        return(out)
      }
    }
  }
  return(rep(NA, (length(index) + 6)*3))
}

prepDat <- function(cor.stats, spec.dat){
  dats <- do.call(rbind, cor.stats)
  out <- data.frame(dats)
  out$Site <- sapply(strsplit(names(cor.stats), "\\."),
                     function(x) x[1])
  out$Date <-  sapply(strsplit(names(cor.stats), "\\."),
                      function(x) x[2])
  out$Year <- format(as.Date(out$Date), "%Y")
  rownames(out) <- NULL
  return(out)
}

vaznull.fast <- function(web) {
  rs.p <- rowSums(web)/sum(web)
  cs.p <- colSums(web)/sum(web)
  P <- P1 <- tcrossprod(rs.p, cs.p)
  finalmat <- matrix(0, nrow(web), ncol(web))
  n.int.finalmat <- 0
  while (n.int.finalmat < sum(dim(web))) {
    sel <- sample(1:length(web), 1, prob = P, replace = TRUE)
    selc <- floor((sel - 1)/(dim(web)[1])) + 1
    selr <- ((sel - 1)%%dim(web)[1]) + 1
    if (sum(finalmat[, selc]) == 0 | sum(finalmat[selr,
                                                  ]) == 0) {
      finalmat[sel] <- 1
      P[sel] <- 0
    }
    n.int.finalmat <- sum(rowSums(finalmat) > 0) + sum(colSums(finalmat) >
                                                         0)
  }
  conn.remain <- sum(web > 0) - sum(finalmat > 0)
  if (conn.remain > 0) {
    if (length(which(finalmat == 0)) == 1) {
      add <- which(finalmat == 0)
    }
    else {
      add <- sample(which(finalmat == 0), conn.remain,
                    prob = P1[finalmat == 0])
    }
    finalmat[add] <- 1
  }
  int.remain <- sum(web) - sum(finalmat)
  if (int.remain > 0) {
    add <- sample(which(finalmat > 0),
                  int.remain, prob = P1[which(finalmat >
                                                0)], replace = TRUE)
    finalmat[as.numeric(names(table(add)))] <-
      finalmat[as.numeric(names(table(add)))] +
      table(add)
  }
  finalmat
}

makeComp <- function(data, metrics, comparison="log") {
  ## takes the sex-level network traits of the pollinators and calculates
  ## the log ratio between males and females within each species and 
  ## network. Within each network, it also drops any species where only
  ## one sex was caught.
  
  som <- mclapply(data, function(x){
    
    #narrow down to each unique Site+Year combo
    col.ls <- lapply(unique(x$SiteYr),function(y){
      sites <- filter(x, x$SiteYr == y)
      
      #narrow down to each unique species within each SiteYr
      spp <- lapply(unique(sites$Sp), function(z){
        sp <- filter(sites,sites$Sp == z)
        
        #If there are two sexes of that species, get the log ratios of each 
        #metric specified in the metrics argument vector     
        if (length(sp$sex)==2) {
          ratios <- lapply(metrics, function(a){
            if(comparison == "log"){
              ratio <- log10(sp[,a][sp$sex == "m"]
                             /sp[,a][sp$sex == "f"])
            } else {
              diff <- sp[,a][sp$sex == "m"] - sp[,a][sp$sex == "f"]
            }
            #  }
          })
          ratios <- data.frame(ratios)
          colnames(ratios) <- metrics
          ratios$species <- z
          ratios$SiteYr <- y
        } 
        else {
          ratios <- NULL
        }
        return(ratios)
      })
      spp <- spp[!sapply(spp, is.null)]
      return(do.call(rbind, spp))
    })
    return(do.call(rbind, col.ls))
  },mc.cores=cores)
  return(som)
}


calcNullProp <- function(data, metrics, zscore=TRUE) {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    data[[x]]$SpSiteYr <- paste(data[[x]]$species,
                                data[[x]]$SiteYr,
                                sep="_")
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(dist.df$SpSiteYr), function(y) {
    sp <- filter(dist.df, dist.df$SpSiteYr == y)
    obs <- filter(sp, sp$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      if (zscore == TRUE){
        metZ <- scale(sp[,z],center = TRUE, scale = TRUE)
        #metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1]) 
        #gotta think about NA treatment here too. I kinda think its
        #getting rid of stuff again. Though this may be fixed when
        #I fix stuff above
      }else{
        metprop <- sum(sp[,z] <= obs[,z]) / length(sp$species)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SpSiteYr <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}

calcNullProp50 <- function(data, metrics, zscore=TRUE) {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    data[[x]]$SpSiteYr <- paste(data[[x]]$species,
                                data[[x]]$SiteYr,
                                sep="_")
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(dist.df$SpSiteYr), function(y) {
    sp <- filter(dist.df, dist.df$SpSiteYr == y)
    obs <- filter(sp, sp$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      if (zscore == TRUE){
        metZ <- scale(sp[,z],center = TRUE, scale = TRUE)
        #metZobs <- ifelse(is.nan(metZ[1]),0,metZ[1]) 
        #gotta think about NA treatment here too. I kinda think its
        #getting rid of stuff again. Though this may be fixed when
        #I fix stuff above
      }else{
        metprop <- (sum(sp[,z] < obs[,z]) + sum(sp[,z] == obs[,z])/2)  / length(sp$species)
      }
    })
    mets <- data.frame(mets)
    colnames(mets) <- metrics
    mets$SpSiteYr <- y
    return(mets)
  },mc.cores=cores)
  
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}


overallTest <- function(prop.dist, metrics, tails = 1, zscore = TRUE) {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2). Can alternatively
  ## give the proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time. 
  
  alpha <- lapply(metrics, function(x){
    clean <- !is.na(prop.dist[,x])
    clean.df <- prop.dist[clean,]
    if(tails == 1){
      if(zscore == TRUE) {
        sum(1.645 <= clean.df[,x]) / length(clean.df[,x])
      } else{
        sum(0.05 >= clean.df[,x]) / length(clean.df[,x])
      }
    } else {
      if(zscore == TRUE) {
        sum(1.96<=clean.df[,x] | -1.96>=clean.df[,x]) / length(clean.df[,x])
      } else{
        sum(0.025>=clean.df[,x] | 0.975<=clean.df) / length(clean.df[,x])
      }
    }
  })
  alpha <- data.frame(alpha)
  colnames(alpha) <- metrics
  return(alpha)
}

spLevelTest <- function(prop.dist, metrics,zscore=TRUE, tails=1) {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2) for each
  ## species in the dataset separately. Can alternatively give the 
  ## proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time
  
  prop.dist$Sp <- gsub( "_.*$", "", prop.dist$SpSiteYr)
  spp <- lapply (unique(prop.dist$Sp),function(x){
    sp <- filter(prop.dist, prop.dist$Sp == x)
    sp.sig <- overallTest(sp, metrics, zscore=zscore, tails=tails)
    sp.sig$Sp <- x
    return(sp.sig)
  })
  spp.sig <- do.call(rbind,spp)
  return(spp.sig)
}



genNullDist <- function(data, metrics, mean.by="SpSiteYr",zscore=TRUE) {
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    data[[x]]$sim <- rep.int(x,times=length(data[[x]]$SiteYr))
    data[[x]]$SpSiteYr <- paste0(data[[x]]$species,data[[x]]$SiteYr)
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #extract an average value w/in each [mean.by] (e.g., mean.by="sim" 
  #would calculate a single mean value w/in each simulation across sites
  dist.build <- mclapply(unique(dist.df[,mean.by]), function(y) {
    net <- filter(dist.df, dist.df[,mean.by] == y)
    obs <- filter(net, net$sim == 1)
    
    #calculate the proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      if (zscore == TRUE){
        z.dist <- scale(net[,z],center=T,scale=T)
        exp.obs <- z.dist[1,]
      }else{
        mean.value <- mean(net[,z],na.rm=T)
        exp.obs <- mean.value - obs[,z]
      }
      return(exp.obs)
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