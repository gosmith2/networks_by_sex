## Functions used during analyses of the Network by Sex project



fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)
  
  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))
  
  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}

dat.clean <- function(spec.dat) {
  spec.dat$GenusSpecies <- fix.white.space(paste(spec.dat$Genus,
                                                 spec.dat$Species,
                                                 spec.dat$SubSpecies))
  
  spec.dat$PlantGenusSpecies <-  fix.white.space(paste(spec.dat$PlantGenus,
                                                       spec.dat$PlantSpecies,
                                                       spec.dat$PlantVar,
                                                       spec.dat$PlantSubSpecies))
  
  spec.dat$Int <-  fix.white.space(paste(spec.dat$GenusSpecies,
                                         spec.dat$PlantGenusSpecies))
  spec.dat$IntGen <-  fix.white.space(paste(spec.dat$Genus,
                                            spec.dat$PlantGenus))
  return(spec.dat)
}

dat.dates <- function(spec.dat) {
  spec.dat$Date <- as.Date(spec.dat$Date, format='%m/%d/%y')
  spec.dat$Doy <- as.numeric(strftime(spec.dat$Date, format='%j'))
  spec.dat$Year <- as.numeric(format(spec.dat$Date,'%Y'))
  return(spec.dat)
}

dat.rm.blanks <- function(spec.dat) {
  spec.dat <- spec.dat[spec.dat$PlantGenusSpecies != "",]
  spec.dat <- spec.dat[spec.dat$GenusSpecies != "",]
  return(spec.dat)
}

samp2site.spp <- function(site, spp, abund, FUN=mean) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
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



## calculates the species roles from a network and returns a dataframe
## with site status and ypr
## takes networks and specimen data

calcSpec <- function(nets, indiv = 1, lvl = "SpSiteYr",index=c('degree','d','weighted.closeness')){
  ## applies specieslevel from bipartite to networks
  species.lev <- lapply(nets, function(x){
    #calculate values
    sl <- specieslevel(x,index=index)
    
    #keep only species with observations >= indiv
    x <- as.data.frame(x)
    sums <- lapply(x, sum)
    cullNames <- names(sums[sums >= indiv])

    sl$`higher level` <- sl$`higher level`[cullNames,]
    
    return(sl)
  })
  #sum up under species (rowsums or such)
  
  ## extract the values and make a dataframe
  if(lvl=="SpSiteYr") {
    specs <- mapply(function(a, b)
      getSpec(species.lev = a,
              names.net = b,
              seps="[.]"),
      a = species.lev,
      b = names(nets),
      SIMPLIFY = FALSE)
  } else {
    specs <- mapply(function(a, b)
      getSpecSpYr(species.lev = a,
              names.net = b,
              seps="[.]"),
      
      SIMPLIFY = FALSE)
  }
  
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

calcNets <- function(nets, metrics,null.sim = F,null.num = 100,null.fun = vaznull.fast){
  ## calculates netowrk-level metrics from networks
  
  #internal null simulaitons off - returns raw metric values
  if(null.sim==F){
    net.lev <- lapply(nets, function(x){
      #calculate values
      sl <- networklevel(x, index=metrics)
      return(sl)
    })
    
  #internal null simulations on - returns zscores of metrics relative to null networks with random interactions
  } else{
    net.lev <- lapply(nets,function(x){
      sl <- calcNetworkMetrics(x,null.num,index=metrics,null.method=null.fun)
      sl <- sl[grep('z',names(sl))]
    })
  }
  
  
  ## extract the values and make a dataframe
  stats <- as.data.frame(do.call(rbind, net.lev))
  stats$SiteYr <- rownames(stats)
  rownames(stats) <- NULL
  return(stats)
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


plant.shuffler <- function(spec.data, trt){
  col.ls <- lapply(unique(spec.data$SiteYr),function(x){
    net <- filter(spec.data,SiteYr==x)
    sp.col <- lapply(unique(net$GenusSpecies), function(y){
      sp<-filter(net,net$GenusSpecies==y)
      
      sp <- filter(sp,sp$Sex %in% c("f","m"))
      
      if(length(unique(sp$Sex)) !=1){
        ##randomize the sexes in a new column, paste together with species name
        #sp$MixSex<-sample(sp$Sex,replace=FALSE)
        
        ##replace plant column with one where males and females share same visitation vector.
          ##this method is random; can't think of a great way right now to do a more concious 
          ##50/50 split that wouldn't run into issues with plants visited only once
          ##This should definitely smooth things out to some degree though
        if(trt=='same'){
          SamePlant <- sample(sp$PlantGenusSpecies, replace=FALSE)
          sp$PlantGenusSpecies <- SamePlant
          
        #replace plant column with one where visitation has no overlaps between males and females
        } else if(trt=='diff'){
  
          if(any(unique(sp$PlantGenusSpecies[sp$Sex=='m']) %in% unique(sp$PlantGenusSpecies[sp$Sex=='f'])) 
             & length(unique(sp$PlantGenusSpecies))>1){
            
            ##make list of plants visited by the species
            plants <- unique(sp$PlantGenusSpecies)
            
            over <- unique(sp$PlantGenusSpecies[sp$Sex=='f'])[unique(sp$PlantGenusSpecies[sp$Sex=='f']) 
                                                              %in% unique(sp$PlantGenusSpecies[sp$Sex=='m'])]
            
            #solve issue where all plant species overlap
            if(length(over)<length(unique(sp$PlantGenusSpecies))){
              
              ##So that things don't break, first assign the overlapping plants to the sex with lower degree
              ##lower than and equal biases slightly towards females having higher degree
              if(length(unique(sp$PlantGenusSpecies[sp$Sex=='m']))
                 >=length(unique(sp$PlantGenusSpecies[sp$Sex=='f']))){
                
                fPlants <- unique(sp$PlantGenusSpecies[sp$Sex=='f'])
                mPlants <- plants[!plants %in% fPlants]
                
              } else{
                mPlants <- unique(sp$PlantGenusSpecies[sp$Sex=='m'])
                fPlants <- plants[!plants %in% mPlants]
              }
              
            } else {
              fPlants <- sample(unique(sp$PlantGenusSpecies[sp$Sex=='f']),size = length(unique(sp$PlantGenusSpecies))/2)
              mPlants <- plants[!plants %in% fPlants]
              
            }
            ##assign individuals to those plants
            #sp$DiffPlant <- sp$PlantGenusSpecies
            sp$PlantGenusSpecies[sp$Sex=='m' & sp$PlantGenusSpecies %notin% mPlants] <- sample(mPlants,size=length(sp$PlantGenusSpecies[sp$Sex=='m' & sp$PlantGenusSpecies %notin% mPlants]),replace=T)
            sp$PlantGenusSpecies[sp$Sex=='f' & sp$PlantGenusSpecies %notin% fPlants] <- sample(fPlants,size=length(sp$PlantGenusSpecies[sp$Sex=='f' & sp$PlantGenusSpecies %notin% fPlants]),replace=T)
            
          } else {
            sp$PlantGenusSpecies <- sp$PlantGenusSpecies
          }
        } else {
          sp$PlantGenusSpecies <- sp$PlantGenusSpecies
        }

      } else {
        sp$PlantGenusSpecies <- sp$PlantGenusSpecies
      }
      
      return(sp)
      })
    })
  return(do.call(rbind,unlist(col.ls,recursive=FALSE)))
}


`%notin%` <- Negate(`%in%`)

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

ran.gen<-function(spec.data,iterations,cores, level="SpSiteYr",boot=FALSE,bootnum=c(2,3,1)){
  #generates simulated datasets by randomizing males and females within species
  
  #setup: add column to actual observation df
  spec.data$MixSex<-spec.data$Sex
  spec.data$GenusSpeciesMix<-paste(spec.data$GenusSpecies,
                                   spec.data$MixSex,
                                   sep="_")
  
  #Keep only networks where at least one sp has both sexes
  spec.keeps <- removeNets(spec.data)
  
  
  #optional chunk: subsample each sp:sex combo w/in each network
  if(boot==TRUE){
    data_boot <- lapply(c(1:3),function(a){
      data_sites <- subset(spec.keeps,spec.keeps$dataset==unique(spec.keeps$dataset)[[a]])
      
      sites_boot <- mclapply(unique(data_sites$SiteYr),function(x){
        site <- subset(data_sites,data_sites$SiteYr==x)
        
        booted <- lapply(unique(site$GenusSpeciesSex),function(y){
          sp <- subset(site,site$GenusSpeciesSex==y)
          len <- length(sp$GenusSpeciesSex)
          bootsamp <- sample_n(sp,(len*bootnum[[a]]),replace=T)
          return(bootsamp)
          })
      
        booted <- do.call(rbind,booted)
        return(booted)
        })
      return(do.call(rbind,sites_boot))
    })
    spec.keeps <- do.call(rbind,data_boot)
  }
  
  #innitiate list
  destList<-list() 
  destList[[1]]<-spec.keeps
  
  #scramble sex column within each species within each network
  it.vec<-1:iterations
  randoms<-mclapply(it.vec, function(z){
    if(level=="SpSiteYr") {
      spec.keeps<-ran.sex(spec.keeps)
      spec.keeps$GenusSpeciesMix<-paste(spec.keeps$GenusSpecies,
                                        spec.keeps$MixSex,
                                        sep="_")
    } else {
      spec.keeps<-ran.sex.SpYr(spec.keeps)
      spec.keeps$GenusSpeciesMix<-paste(spec.keeps$GenusSpecies,
                                        spec.keeps$MixSex,
                                        sep="_")
    }
    print(z)
    return(spec.keeps)
  }, mc.cores = cores)
  destList<-c(destList,randoms)
  #check that the sexes actually mixed
  
  ifelse(any(destList[[2]]$Sex!=destList[[2]]$MixSex),
         print("SUCCESS: the sexes mixed correctly"),
         print("WARNING: the sexes did not mix correctly"))
  
  return(destList)
}





## function to simulate 1 null, and calculate statistics on it
calcNullStat <- function(dat.web,
                         null.fun,...) {
  sim.web <- null.fun(dat.web)
  return(calcMetric(sim.web,...))
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
calcNetworkMetrics <- function (dat.web, N, null.method,
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
                                            null.fun= null.method,
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

basnull <- function (web) 
  #modified / simplified from bipartite::shuffle.web, which itself is based on Fortuna&Bascompte 2006 null model 1 
{
  web <- as.matrix(web)
  web <- empty(web)
  if (dim(web)[1] > dim(web)[2]) {
    long = TRUE
    web = t(web)
  } else {
    long = FALSE
  }
  shuffle <- function(web) {
    dimdiff <- dim(web)[1] - dim(web)[2]
    if ((length(diag(web)) + dimdiff) > sum(as.vector(web) > 
                                            0)) 
      stop("Too few entries in the web: less interactions than length of web diagonal.")
    out <- web
    out[, ] <- 0
    shuf <- sample(as.vector(web))
    nozero.index <- which(shuf != 0)
    diag(out) <- shuf[nozero.index[1:length(diag(out))]]
    xdiag.index = matrix(ncol = 2, nrow = max(dim(web)) - 
                             min(dim(web)))
    xdiag.index[, 1] = sample(c(1:(min(dim(web)))), size = max(dim(web)) - 
                                  min(dim(web)), replace = TRUE)
    xdiag.index[, 2] = (min(dim(web)) + 1):max(dim(web))
    out[xdiag.index] = shuf[nozero.index[(length(diag(web)) + 
                                              1):(length(diag(web)) + nrow(xdiag.index))]]
    gone <- sum(out > 0)
    remains <- shuf[nozero.index[-c(1:gone)]]
    option <- which(out == 0, arr.ind = TRUE)
    out[option[sample(dim(option)[1], length(remains)), 
                 ]] <- remains
    colnames(out) <- rownames(out) <- NULL
    if (long) {
        out = t(out)
      }
      out
  }
  shuffle(web)
}


makeComp <- function(data, metrics, comparison="log") {
  ## takes the sex-level network traits of the pollinators and calculates
  ## the log ratio (or difference) between males and females within each species and 
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

calcNullProp50 <- function(data, metrics, zscore=TRUE,level ="species") {
  ## calculates the zscore of the observed difference between male and
  ## female pollinators of each species within the larger distribution
  ## of that difference across all iterations. Can alternatively give
  ## the proportion of iterations greater than the observed value
  
  #give each element of the long list a unique number
  sim.vec <- seq(1:length(data))
  named.ls <- mclapply(sim.vec, function(x) {
    data[[x]]$sim <- rep.int(x,times=length(data[[x]][,1]))
    if(level=="species"){
      data[[x]]$SpSiteYr <- paste(data[[x]]$species,
                                  data[[x]]$SiteYr,
                                  sep="_")
    }
      
    return(data[[x]])
  },mc.cores=cores)
  
  #combine all the simulation iterations together into single element
  dist.df <- do.call(rbind, named.ls)
  
  #combine values for sp+yr+site
  if(level == "species") {
    col <- "SpSiteYr"
  } else if(level=='combo'){
    col <- 'combo'
  } else {
    col <- "SiteYr"
  }
  sigLevel <- mclapply(unique(dist.df[,col]), function(y) {
    sp <- filter(dist.df, dist.df[,col] == y)
    obs <- filter(sp, sp$sim == 1)
    
    print(1)
    #calculate the zscore of the observed difference, 
    #or proportion of simulations <= observed
    mets <- lapply(metrics, function(z) {
      if (zscore == TRUE){
        metZ <- (sp[1,z] - mean(sp[,z]))/
                 (sd(sp[,z])+10^-10)
      }else{
        metprop <- (sum(sp[,z] < obs[,z]) + sum(sp[,z] == obs[,z])/2)  / length(sp[,z])
      }
    })
    mets <- data.frame(t(unlist(mets)))
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
        y <- sum(1.645 <= clean.df[,x]) / length(clean.df[,x])
      } else{
        y <- sum(0.05 >= clean.df[,x]) / length(clean.df[,x])
      }
    } else if(tails == 2){
      if(zscore == TRUE) {
        y <- sum(1.96<=clean.df[,x] | -1.96>=clean.df[,x]) / length(clean.df[,x])
      } else{
        y <- sum(0.025>=clean.df[,x] | 0.975<=clean.df[,x]) / length(clean.df[,x])
      }
    } else {
      if(zscore == TRUE) {
        y <- sum(-1.645 >= clean.df[,x]) / length(clean.df[,x])
      } else{
        y <- sum(0.95 <= clean.df[,x]) / length(clean.df[,x])
      }
    }
  })
  alpha <- data.frame(alpha)
  colnames(alpha) <- metrics
  return(alpha)
}

spLevelTest <- function(prop.dist, metrics,zscore=TRUE, tails=1,level="GenusSpecies") {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2) for each
  ## species in the dataset separately. Can alternatively give the 
  ## proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time
  spp <- lapply (unique(prop.dist[,level]),function(x){
    sp <- filter(prop.dist, prop.dist[,level] == x)
    sp.sig <- overallTest(sp, metrics, zscore=zscore, tails=tails)
    sp.sig$level <- x
    return(sp.sig)
  })
  spp.sig <- do.call(rbind,spp)
  return(spp.sig)
}
#May not use; not in core code



distComp <- function(dataset, compMethod, indiv = 1) {
  #narrow down to each randomized simulation
  sim.vec <- seq(1:length(dataset))
  simDistance <- mclapply(sim.vec, function(x) {
    sim <- dataset[[x]]
    
    #narrow down to each network in that simulation
    net.vec <- seq(1:length(sim))
    netDist <- lapply(net.vec, function(y) {
      net <- sim[y]
      
      #Make list of pollinator nodes in that network
      nodeList <- colnames(net[[1]])
      
      #remove nodes where sex wasn't defined
      nodeListClean <- unlist(lapply(nodeList,function(a){
        if(gsub( "^.*_", "", a) == "m" |gsub( "^.*_", "", a) == "f"){
          return(a)
        }
      }))
      
      #generate list of only species where both sexes are present
      spList <- gsub( "_.*$", "", nodeListClean)
      sexList <- unlist(spList[duplicated(spList)])
      
      #generate list of only species where both sexes have a minimum number of individuals
      minList <- lapply(sexList, function(b){
        min <- min(colSums(net[[1]][,c(paste0(b,"_m"),paste0(b,"_f"))]))
      })
      minDf <- data.frame(Sp = sexList,min = unlist(minList)) 
      targetList <- minDf$Sp[minDf$min>=indiv]
      
      #for the target species, calculate the distance index
      if(length(targetList)==0){
        netVal <- data.frame(GenusSpecies=NA,distance=NA,SiteYr=names(net))
        
      } else {
        spDist <- lapply(targetList, function(z){
          
          comp <- net[[1]][,c(paste0(z,"_m"),paste0(z,"_f"))]
          dist <- as.matrix(vegdist(t(comp),method=compMethod))
          data.frame(GenusSpecies=z,distance=dist[2,1])
        })
        
        netVal <- do.call(rbind,spDist)
        netVal$SiteYr <- rep(names(net))
        return(netVal)
      }
    })
    simDist <- do.call(rbind, netDist)
    simDist$sim <- rep(x)
    return(simDist)
  },mc.cores=cores)
  allDist <- do.call(rbind, simDistance)
  return(allDist)
}

distTest <- function(zDist, tails = 1) {
  ## calculates the proportion of z scores above a threshold based
  ## on tails (5% for 1 tail, 2.5% each direction for 2). Can alternatively
  ## give the proportion of sp+site+yr observations whose iterations
  ## differed from the observed over 95% of the time. 
  
  if(tails == 1){
    y <- sum(1.645 <= zDist$distance) / length(zDist$distance)
    
  } else {
    y <- sum(1.96<=zDist$distance | -1.96>=zDist$distance) / length(zDist$distance)
  }
  return(y)
}

calcDistZ <- function(data, level, zscore = TRUE) {
  #combine values for sp+yr+site
  sigLevel <- mclapply(unique(data[[level]]), function(x) {
    lev <- filter(data, data[[level]] == x)
    obs <- filter(lev, lev$sim == 1)

    #calculate the zscore of the observed difference, 
    #or proportion of simulations <= observed
    mets <- 
      if (zscore == TRUE){
        metZ <- (lev[1,'distance'] - mean(lev[,'distance']))/
          (sd(lev[,'distance'])+10^-10)
      }else{
        metprop <- (sum(lev[,'distance'] < obs[,'distance']) + sum(lev[,'distance'] == obs[,'distance'])/2)  / length(lev[,level])
      }
    mets <- data.frame(distanceZ=mets,Level=x)
    return(mets)
  })
  
  #bind these all together
  sig.dist <- do.call(rbind,sigLevel)
  return(sig.dist)
}


traitFrame <- function(diffs, dist,level='SpSiteYr') {
  #extract species from SpSiteYr
  if(level=='SpSiteYr'){
    dist %>%
      mutate(obs = gsub("\\.", "_", Level)) ->
      dist1
    
    #add in distance index proportions, matching by SpSiteYr
    diffs$distance <- dist1$distanceZ[match(diffs$SpSiteYr,
                                            dist1$obs)]
    diffs$Site <- word(diffs$SpSiteYr,2,3,sep="_")
    return(diffs)
    
  } else {
    diffs %>%
      mutate(obs=combo) ->
      diff1
    dist %>%
      mutate(obs=combo) ->
      dist1
    
    diff1$distance <- dist1$distanceZ[match(diff1$obs,
                                            dist1$obs)]
    return(diff1)
  }
  
  return(diff1)
}



zscore_avger <- function(data,metrics){
  data$year <- word(data$SpSiteYr,3,sep='_')
  data$species <- word(data$SpSiteYr,1,sep='_')
  data$Site <- word(data$SpSiteYr,2,sep='_')
  data$SiteYr <- paste(data$Site,data$year,sep="_")
  data$dataset <- spec_all$dataset[match(data$SiteYr,sub(" ","_",spec_all$SiteYr))]
  data$combo <- paste(data$species,data$year,data$dataset,sep="_")
  
  avg <- lapply(unique(data$combo),function(x){
    combo <- data[data$combo==x,]
    
    mets <- do.call(cbind,lapply(metrics, function(y){
      metavg <- mean(combo[,y],na.rm=T)
      metavg <- data.frame('met' = metavg) 
      names(metavg) <- y
      return(metavg)
    }))
    
    mets$combo <- x
    
    return(mets)
  })
  
  avgs <- do.call(rbind,avg)
  
  return(avgs)
}

distZ_avger <- function(data){
  data$year <- substr(data$Level,str_length(data$Level)-3,str_length(data$Level))
  data$species <- gsub("_.*$", "", data$Level)
  data$SiteYr <- sub('.*_',"",data$Level)
  data$dataset <- spec_all$dataset[match(data$SiteYr,sub(" ",".",spec_all$SiteYr))]
  data$combo <- paste(data$species,data$year,data$dataset,sep="_")
  
  avg <- lapply(unique(data$combo),function(x){
    combo <- data[data$combo==x,]
    
    dist <- data.frame('distanceZ'=mean(combo$distanceZ))
    
    dist$combo <- x
    
    return(dist)
  })
  
  avgs <- do.call(rbind,avg)
  
  return(avgs)
  
}

t.tester <- function(data,metric){
  
  mets <- lapply(metric,function(x){
    t <- t.test(data[,x])
    df <- data.frame('metric' = x, 'tval' = round(t$statistic[[1]],digits=4), 'pval' = round(t$p.value,digits=4),'mean' = mean(data[,x],na.rm=T))
    return(df)
  })
  
  do.call(rbind,mets)
}

lme.maker <- function(data,metric,random.eff){
  mets <- lapply(metric,function(x){
    browser()
    lme.met <- lme(data[[,x]]~1,random=~1|data[,random.eff])
    lmes <- summary(lme.met)
    df <- as.data.frame(lmes$tTable)
    df$metric <- x
    df$mean <- mean(data[,x],na.rm=T)
    return(df)
  })
  do.call(rbind,mets)
}

wholeNet <- function(sampData){
  obs <- sampData[sampData$Sex=="f",]
  agg.spec <- aggregate(list(abund=obs$GenusSpecies),
                        list(GenusSpecies=
                               obs$GenusSpecies,
                             PlantGenusSpecies=
                               obs$PlantGenusSpecies),
                        length)
  net <- samp2site.spp(agg.spec$PlantGenusSpecies, 
                       agg.spec$GenusSpecies,
                       agg.spec$abund)
  return(net)
}

calcRareDeg <- function(sampData, name){
  net <- wholeNet(sampData)
  rare <- apply(net,2,chao1)
  rare.df <- data.frame(GenusSpecies = names(rare),
                        r.deg = rare,
                        row.names = NULL)
  names(rare.df) <- c("GenusSpecies",name)
  return(rare.df)
}

zscore_rare_joiner <- function(nodez,distz,networks,meta,metric, model){
  nodez$dist <- distz$distanceZ[match(nodez$SpSiteYr,gsub("\\.","_",distz$Level))]
  nodez$GenusSpecies <- word(nodez$SpSiteYr,1,sep = "_")
  nodez$Genus <- word(nodez$GenusSpecies,1,sep=" ")
  nodez$SiteYr <- word(nodez$SpSiteYr,2,3,sep='_')
  nodez$diffs <- unlist(lapply(nodez$SpSiteYr,function(x){
    sp <- nodez[nodez$SpSiteYr==x,]
    val <- sum(sp[metric]>1.96|sp[metric]<(-1.96), na.rm=T)
    return(val)
  })
  )
  nodez$rareDeg <- rare.df$r.deg1[match(nodez$GenusSpecies,rare.df$GenusSpecies)]
  nodez$trials <- 4
  nodez$plants <- unlist(lapply(nodez$SiteYr,function(x){
    nets <- networks[[1]]
    net <- nets[[gsub("_",".",x)]]
    dim <-dim(net)
    return(dim[1])
  }))
  nodez$dataset <- as.factor(meta$dataset[match(nodez$SiteYr,gsub(" ","_",meta$SiteYr))])
  nodez$Family <- meta$Family[match(nodez$GenusSpecies,meta$GenusSpecies)]
  nodez$model <- as.factor(rep(model))
  return(nodez)
}

zscore_net_joiner <- function(netz,networks,meta, model){
  netz$plants <- unlist(lapply(netz$SpSiteYr,function(x){
    nets <- networks[[1]]
    net <- nets[[x]]
    dim <-dim(net)
    return(dim[1])
  }))
  netz$dataset <- as.factor(meta$dataset[match(netz$SpSiteYr,gsub(" ",".",meta$SiteYr))])
  netz$model <- as.factor(rep(model))
  return(netz)
}

  

calcMetric <- function(dat.web, ...) {
  ## calculates modularity
  dat.web <- as.matrix(empty(dat.web))
  ## matrix of all the same number
  if(min(dat.web) == max(dat.web)){
    return(c(mets=rep(0, 1),
             mod.met=rep(0,2)))
  } else{
    mets <-  networklevel(dat.web, ...)
  }
  
  ## the functional redundancy function takes a matrix of sites and
  ## species, and a trait matrix whwere the rownames of the traits
  ## match the column names of the site by species matric. In our
  ## case, the plants and pollinators are the "species"
  ## respectively, and their traits are their interaction partners.
  
  ## create names for later use in site x species matrix
  rownames(dat.web) <- 1:nrow(dat.web)
  colnames(dat.web) <- 1:ncol(dat.web)
  
  ## site by species matrix where there is only one "site" since the
  ## networks are site specific, and the columns are the
  ## species.
  
  ## abundance weighted site by species matrices
  plants <- matrix(rowSums(dat.web),  nrow=1)
  pols <- matrix(colSums(dat.web),  nrow=1)
  
  colnames(plants) <- rownames(dat.web)
  colnames(pols) <- colnames(dat.web)
  
  ## pull out Functional redundancy score based on: de Bello, F.;
  ## Leps, J.; Lavorel, S. & Moretti, M. (2007). Importance of
  ## species abundance for assessment of trait composition: an
  ## example based on pollinator communities. Community Ecology, 8,
  ## 163:170. and functional complementarity score based on Rao,
  ## C.R. (1982). Diversity and dissimilarity coefficients: a
  ## unified approach. Theoretical Population Biology, 21, 24:43.
  
  return(mets)
}



extinction1 <- function (web, participant = "both", method = "degree", 
          ext.row = NULL, ext.col = NULL) {
  partis <- c("lower", "higher", "both")
  partis.match <- pmatch(participant, partis)
  if (is.na(partis.match)) 
    stop("Choose participant: lower/higher/both.\n")
  meths <- c("random", "abundance", "degree", 
             "external")
  meths.match <- pmatch(method, meths)
  if (is.na(meths.match)) 
    stop("Choose extinction method: random/abundance/degree.\n")
  nr <- NROW(web)
  nc <- NCOL(web)
  if (partis.match == 3 & meths.match == 1) 
    partis.match <- sample(2, 1)
  if (meths.match == 1) {
    rexcl <- sample(nr, 1)
    cexcl <- sample(nc, 1)
    if (partis.match == 1) 
      web[rexcl, ] <- 0
    if (partis.match == 2) 
      web[, cexcl] <- 0
  }
  if (meths.match == 2) {
    web <- web[sample(1:nrow(web)), sample(1:ncol(web)), 
               drop = FALSE]
    rseq <- order(rowSums(web))
    cseq <- order(colSums(web))
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) {
      if (word(colnames(web)[cseq[1]],2,sep="_")=="m"){
        web[, cseq[1]] <- 0
      }
      if (word(colnames(web)[cseq[1]],2,sep="_")=='f'){
        sp <- word(colnames(web)[cseq[1]],1,sep="_")
        web[,grepl(sp, colnames(web))]<-0
      }
    }
    if (partis.match == 3) {
      if (min(rowSums(web)) < min(colSums(web))) {
        web[rseq[1], ] <- 0
      }
      else {
        if (min(rowSums(web)) > min(colSums(web))) {
          if (word(colnames(web)[cseq[1]],2,sep="_")=="m"){
            web[, cseq[1]] <- 0
          }
          if (word(colnames(web)[cseq[1]],2,sep="_")=='f'){
            sp <- word(colnames(web)[cseq[1]],1,sep="_")
            web[,grepl(sp, colnames(web))]<-0
          } else {
            web[, cseq[1]] <- 0
          }
        }
        else {
          if (sample(2, 1) == 1) {
            web[rseq[1], ] <- 0
          }
          else {
            if (word(colnames(web)[cseq[1]],2,sep="_")=="m"){
              web[, cseq[1]] <- 0
            }
            if (word(colnames(web)[cseq[1]],2,sep="_")=='f'){
              sp <- word(colnames(web)[cseq[1]],1,sep="_")
              web[,grepl(sp, colnames(web))]<-0
            }
          }
        }
      }
    }
  }
  if (meths.match == 3) {
    if (partis.match == 1) {
      sequ <- rowSums(web > 0)
      which.ex <- which(sequ == max(sequ))
      if (length(which.ex) > 1) {
        ex <- sample(which.ex, size = 1)
      }
      else {
        ex <- which.ex
      }
      web[ex, ] <- 0
    }
    if (partis.match == 2) {
      sequ <- colSums(web > 0)
      which.ex <- which(sequ == max(sequ))
      if (length(which.ex) > 1) 
        ex <- sample(which.ex, size = 1)
      else ex <- which.ex
      if (word(names(ex),2,sep='_')=="m"){
        web[, ex] <- 0
      }
      if (word(names(ex),2,sep='_')=='f'){
        sp <- word(names(ex),1,sep='_')
        web[,grepl(sp, colnames(web))]<-0
      }
      else{
        web[, ex] <- 0
      }
    }
    if (partis.match == 3) {
      if (min(rowSums(web>0)) < min(colSums(web>0))) {
        sequ <- rowSums(web > 0)
        which.ex <- which(sequ == max(sequ))
        if (length(which.ex) > 1) {
          ex <- sample(which.ex, size = 1)
        }
        else {
          ex <- which.ex
        }
        web[ex, ] <- 0
      }
      if (min(rowSums(web>0)) > min(colSums(web>0))){
        sequ <- colSums(web > 0)
        which.ex <- which(sequ == max(sequ))
        if (length(which.ex) > 1) 
          ex <- sample(which.ex, size = 1)
        else ex <- which.ex
        if (word(names(ex),2,sep='_')=="m"){
          web[, ex] <- 0
        }
        if (word(names(ex),2,sep='_')=='f'){
          sp <- word(names(ex),1,sep='_')
          web[,grepl(sp, colnames(web))]<-0
        }
        else{
          web[, ex] <- 0
        }
      }
      else {
        if (sample(2, 1) == 1) {
          sequ <- rowSums(web > 0)
          which.ex <- which(sequ == max(sequ))
          if (length(which.ex) > 1) {
            ex <- sample(which.ex, size = 1)
          }
          else {
            ex <- which.ex
          }
          web[ex, ] <- 0
        }
        else {
          sequ <- colSums(web > 0)
          which.ex <- which(sequ == max(sequ))
          if (length(which.ex) > 1) 
            ex <- sample(which.ex, size = 1)
          else ex <- which.ex
          if (word(names(ex),2,sep='_')=="m"){
            web[, ex] <- 0
          }
          if (word(names(ex),2,sep='_')=='f'){
            sp <- word(names(ex),1,sep='_')
            web[,grepl(sp, colnames(web))]<-0
          }
          else{
            web[, ex] <- 0
          }
        }
      }
    }
  }
  if (meths.match == 4) {
    rseq <- ext.row
    cseq <- ext.col
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) 
      web[, cseq[1]] <- 0
  }
  if (any(word(colnames(web),2,sep = '_')=='m')){
    males <- colnames(web)[word(colnames(web),2,sep = '_')=='m']
    spTF <- unlist(lapply(males,function(x){
      sp <- word(x,1,sep='_')
      if(length(colnames(web)[grepl(sp,colnames(web))])==1){
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    }))
    spls <- males[spTF]
    web[,spls] <- 0
  }
  return(web)
}

second.extinct1 <- function (web, participant = "both", method = "degree", 
                             nrep = 100, details = FALSE, ext.row = NULL, ext.col = NULL) 
{
  if (participant == "both" & method == "external") 
    stop("Sorry, that won't work. When you specify the sequence, you have to choose one of the two levels. 'both' won't work.")
  if (!is.null(ext.row) & length(ext.row) != NROW(web)) 
    stop("The length of the external row vector is different from the numbers of rows in the network!")
  if (!is.null(ext.col) & length(ext.col) != NCOL(web)) 
    stop("The length of the external col vector is different from the numbers of cols in the network!")
  if (participant == "higher" & method == "external" & 
      is.null(ext.col)) 
    stop("You need to provide an external sequence of extinction for the higher trophic level!")
  if (participant == "lower" & method == "external" & 
      is.null(ext.row)) 
    stop("You need to provide an external sequence of extinction for the lower trophic level!")
  one.second.extinct <- function(web = web, participant = participant, 
                                 method = method, ext.row = ext.row, ext.col = ext.col) {
    dead <- matrix(nrow = 0, ncol = 3)
    colnames(dead) <- c("no", "ext.lower", "ext.higher")
    m2 <- web
    i <- 1
    repeat {
      n <- extinction1(m2, participant = participant, method = method, 
                      ext.row = ext.row, ext.col = ext.col)
      dead <- rbind(dead, c(i, attributes(m2 <- empty(n, 
                                                      count = TRUE))$empty))
      if (participant == "lower" & NROW(m2) < 2) 
        break
      if (participant == "higher" & NCOL(m2) < 2) 
        break
      if (participant == "both" & min(dim(m2)) < 
          2) 
        break
      if (any(dim(n) == 1)) 
        break
      if (method == "external") {
        ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > 
                                                   ext.col[1]] - 1
        ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > 
                                                   ext.row[1]] - 1
        ext.row <- ext.row[-1]
        ext.col <- ext.col[-1]
      }
      i <- i + 1
    }
    dead2 <- rbind(dead, c(NROW(dead) + 1, NROW(m2), NCOL(m2)))
    if (participant == "lower" & method == "degree") {
      if (length(table(dead[, 2])) > 1) 
        dead2[, 2] <- 1
    }
    if (nrow(dead) + 1 != nrow(dead2)) 
      stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
    if (participant == "lower") 
      supposed.length <- NROW(web)
    if (participant == "higher") 
      supposed.length <- NCOL(web)
    if (participant == "both") 
      supposed.length <- NROW(dead2)
    if (NROW(dead2) != supposed.length) {
      missing <- supposed.length - NROW(dead2)
      addit1 <- (NROW(dead2) + 1):(NROW(dead2) + missing)
      addit2n3 <- rep(0, times = missing)
      dead2 <- rbind(dead2, as.matrix(data.frame(addit1, 
                                                 addit2n3, addit2n3)))
    }
    return(dead2)
  }
  if (is.vector(method)) 
    sequence = method
  if (pmatch(method, c("abundance", "random", "degree", 
                       "external")) %in% c(1, 3, 4)) {
    out <- one.second.extinct(web = web, participant = participant, 
                              method = method, ext.row = ext.row, ext.col = ext.col)
  }
  else {
    o <- replicate(nrep, one.second.extinct(web = web, participant = participant, 
                                            method = method, ext.row = ext.row, ext.col = ext.col), 
                   simplify = FALSE)
    if (details) {
      out <- o
    }
    else {
      lengths <- sapply(o, nrow)
      z <- o[[which.max(lengths)]]
      z[, 2:3] <- 0
      for (k in 1:length(o)) {
        nr <- nrow(o[[k]])
        z[1:nr, ] <- z[1:nr, ] + o[[k]]
        rm(nr)
      }
      out <- z/length(o)
      out[, 1] <- 1:max(lengths)
    }
  }
  class(out) <- "bipartite"
  attr(out, "exterminated") <- c("both", "lower", 
                                 "higher")[pmatch(participant, c("both", "lower", 
                                                                 "higher"))]
  out
}

robustness1 <- function (object) 
{
  if (!is(object, "bipartite")) 
    stop("This function cannot be meaningfully applied to objects of this class.")
  N <- colSums(object)
  if (attr(object, "exterminated") == "lower") {
    y <- -object[, 3]
  } else if (attr(object, 'exterminated')=='higher') {
    y <- -object[, 2]
  } else {
    y <- -object[,c(2,3)]
  }
  if (class(y) == "matrix"){
    lapply(c(1,2),function(z){
      y <- y[,z]
      y <- (sum(y) - cumsum(y))/sum(y)
      x <- (object[, "no"]/max(object[, "no"]))
      ext.curve <- splinefun(x, y)
      ext.area <- integrate(ext.curve, 0, 1)
      return(as.numeric(ext.area[[1]]))
    })
  } else {
    y <- (sum(y) - cumsum(y))/sum(y)
    x <- (object[, "no"]/max(object[, "no"]))
    ext.curve <- splinefun(x, y)
    ext.area <- integrate(ext.curve, 0, 1)
    return(as.numeric(ext.area[[1]]))
    
  }
  
}


empty <- function (web, count = FALSE) 
{
  web[is.na(web)] <- 0
  if (NCOL(web) == 1 | NROW(web) == 1) {
    if (NCOL(web) == 1 & NROW(web) != 1) {
      nr <- sum(web > 0)
      nc <- 1
    }
    if (NROW(web) == 1 & NCOL(web) != 1) {
      nc <- sum(web > 0)
      nr <- 1
    }
    if (NROW(web) == 1 & NCOL(web) == 1) {
      nr <- 1
      nc <- 1
    }
  }
  cempty <- which(colSums(web) == 0)
  rempty <- which(rowSums(web) == 0)
  cind <- if (length(cempty) == 0) 
    1:NCOL(web)
  else (1:NCOL(web))[-cempty]
  rind <- if (length(rempty) == 0) 
    1:NROW(web)
  else (1:NROW(web))[-rempty]
  out <- web[rind, cind, drop = FALSE]
  if (count) 
    attr(out, "empty") <- c(`empty rows` = length(rempty), 
                            `empty columns` = length(cempty))
  return(out)
}


boot_sumarizer <- function(iter,dataframes,labels,metrics,tails=2,zscore=T){
  
  boot.ls <- lapply(iter,function(x){
    
    ot <- t(overallTest(dataframes[[x]], metrics, tails =tails,zscore=zscore))
    tt <- t.tester(dataframes[[x]], metrics)
    tt$prop <- ot
    tt$label <- rep(labels[x])
    
    return(tt)
  })
  
  boots <- do.call(rbind,boot.ls)
  return(boots)
}

function (web, index = "ALLBUTDD", level = "both", 
          weighted = TRUE, ISAmethod = "Bluethgen", SAmethod = "Bluethgen", 
          extinctmethod = "r", nrep = 100, CCfun = median, dist = "horn", 
          normalise = TRUE, empty.web = TRUE, logbase = "e", 
          intereven = "prod", H2_integer = TRUE, fcweighted = TRUE, 
          fcdist = "euclidean", legacy = FALSE) 
{
  if (empty.web) {
    web <- empty(web)
  }
  web.e <- empty(web)
  if (NROW(web) < 2 | NCOL(web) < 2) 
    warning("Web is really too small to calculate any reasonable index. You will get the values nonetheless, but I wouldn't put any faith in them!")
  allindex <- c("number of species", "connectance", 
                "web asymmetry", "links per species", "number of compartments", 
                "compartment diversity", "cluster coefficient", 
                "degree distribution", "mean number of shared partners", 
                "togetherness", "C score", "V ratio", 
                "discrepancy", "nestedness", "NODF", 
                "weighted nestedness", "ISA", "SA", 
                "extinction slope", "robustness", "niche overlap", 
                "weighted cluster coefficient", "weighted NODF", 
                "partner diversity", "generality", "vulnerability", 
                "linkage density", "weighted connectance", 
                "Fisher alpha", "interaction evenness", "Alatalo interaction evenness", 
                "Shannon diversity", "functional complementarity", 
                "H2")
  index <- unique(index)
  wrong.name <- which(is.na(pmatch(index, c(allindex, "ALL", 
                                            "ALLBUTDD", "info", "quantitative", 
                                            "binary", "topology", "networklevel"))))
  if (length(wrong.name) > 0) 
    stop("You selected an index that is not available: ", 
         paste(index[wrong.name], collapse = ", "))
  if (length(index) == 1 & !all(index %in% allindex)) {
    index <- switch(index, ALL = allindex, ALLBUTDD = allindex[-which(allindex == 
                                                                        "degree distribution")], info = c("number of species", 
                                                                                                          "connectance", "web asymmetry", "links per species", 
                                                                                                          "number of compartments"), quantitative = c("weighted cluster coefficient", 
                                                                                                                                                      "weighted nestedness", "weighted NODF", 
                                                                                                                                                      "functional complementarity", "partner diversity", 
                                                                                                                                                      "effective partners", "H2", "diversity", 
                                                                                                                                                      "linkage density", "weighted connectance", 
                                                                                                                                                      "niche overlap"), binary = c("connectance", 
                                                                                                                                                                                   "links per species", "nestedness", "mean number of partners", 
                                                                                                                                                                                   "cluster coefficient", "C-score", "Fisher alpha"), 
                    topology = c("connectance", "cluster coefficient", 
                                 "degree distribution", "togetherness", 
                                 "nestedness", "NODF"), networklevel = c("connectance", 
                                                                         "web asymmetry", "links per species", 
                                                                         "number of compartments", "compartment diversity", 
                                                                         "cluster coefficient", "nestedness", 
                                                                         "NODF", "weighted NODF", "ISA", 
                                                                         "SA", "linkage density", "Fisher alpha", 
                                                                         "diversity", "interaction evenness", 
                                                                         "Alatalo interaction evenness", "H2"), 
                    stop("Your index is not recognised! Typo? Check help for options!", 
                         call. = FALSE))
  }
  if (legacy == FALSE) {
    out <- list()
    if ("connectance" %in% index) {
      suppressWarnings(out$connectance <- sum(web > 0)/prod(dim(web)))
    }
    if ("web asymmetry" %in% index) 
      out$"web asymmetry" <- (NCOL(web) - NROW(web))/sum(dim(web))
    if ("links per species" %in% index) {
      L <- sum(web > 0)/sum(dim(web))
      out$"links per species" = L
    }
    if (any(c("number of compartments", "compartment diversity") %in% 
            index)) {
      CD <- function(co) {
        if (co$n.compart > 1) {
          no <- NA
          for (i in 1:co$n.compart) {
            comp <- which(abs(co$cweb) == i, arr.ind = TRUE)
            no[i] <- length(unique(comp[, 1])) + length(unique(comp[, 
                                                                    2]))
          }
          no <- no/sum(dim(web))
          CD <- exp(-sum(no * log(no)))
        }
        else {
          CD <- NA
        }
        CD
      }
      comps <- try(compart(web.e), silent = TRUE)
      if (inherits(comps, "try-error")) {
        ncompart <- compdiv <- NA
      }
      else {
        ncompart <- comps$n.compart
        compdiv <- CD(comps)
      }
      if ("number of compartments" %in% index) 
        out$"number of compartments" <- as.integer(ncompart)
      if ("compartment diversity" %in% index) 
        out$"compartment diversity" <- compdiv
    }
    if ("cluster coefficient" %in% index) {
      cluster.coef <- function(web, full = FALSE, FUN = mean) {
        web <- as.matrix(web)
        Ci.high <- colSums((web > 0))/nrow(web)
        Ci.low <- rowSums((web > 0))/ncol(web)
        CC <- FUN(Ci.high)
        if (full) 
          out <- list(`cluster coefficient` = CC, 
                      `CC values higher` = Ci.high, `CC values lower` = Ci.low)
        else out <- c(`cluster coefficient` = CC)
        out
      }
      out$"cluster coefficient" = as.numeric(cluster.coef(web, 
                                                          FUN = CCfun, full = FALSE))
    }
    if ("nestedness" %in% index) {
      nest <- try(nestedtemp(web)$statistic, silent = TRUE)
      out$nestedness <- ifelse(inherits(nest, "try-error"), 
                               NA, nest)
    }
    if ("NODF" %in% index) {
      NODF <- try(unname(nestednodf(web, order = TRUE, 
                                    weighted = FALSE)$statistic[3]), silent = TRUE)
      out$NODF <- if (inherits(NODF, "try-error")) 
        NA
      else NODF
    }
    if ("weighted nestedness" %in% index) {
      wine.res <- try(wine(web.e, nreps = nrep)$wine, silent = TRUE)
      out$"weighted nestedness" <- if (!inherits(wine.res, 
                                                 "try-error")) {
        wine.res
      }
      else {
        NA
      }
    }
    if ("weighted NODF" %in% index) {
      wNODF <- try(unname(nestednodf(web, order = TRUE, 
                                     weighted = TRUE)$statistic[3]), silent = TRUE)
      out$"weighted NODF" <- if (inherits(wNODF, 
                                          "try-error")) 
        NA
      else wNODF
    }
    if (any(c("ISA", "interaction strength asymmetry", 
              "dependence asymmetry") %in% index)) {
      depL <- web.e/matrix(rowSums(web.e), nrow = NROW(web.e), 
                           ncol = NCOL(web.e), byrow = FALSE)
      depH <- web.e/matrix(colSums(web.e), nrow = NROW(web.e), 
                           ncol = NCOL(web.e), byrow = TRUE)
      if (ISAmethod == "Bascompte" & "ISA" %in% 
          index) {
        out$"dependence asymmetry" = mean(abs(depL - 
                                                depH)/pmax(depL, depH), na.rm = TRUE)
      }
      if (ISAmethod == "Bluethgen" & "ISA" %in% 
          index) {
        web2 <- web
        web2[, which(colSums(web) == 1)] <- 0
        web2[which(rowSums(web) == 1), ] <- 0
        rowsummat <- matrix(rowSums(web2), nrow = NROW(web2), 
                            ncol = NCOL(web2), byrow = FALSE)
        colsummat <- matrix(colSums(web2), nrow = NROW(web2), 
                            ncol = NCOL(web2), byrow = TRUE)
        depL <- web2/rowsummat
        depH <- web2/colsummat
        depL[depL <= 0] <- NA
        depH[depH <= 0] <- NA
        depLprime <- (depL - 1/rowsummat)/(1 - 1/rowsummat)
        depHprime <- (depH - 1/colsummat)/(1 - 1/colsummat)
        out$"interaction strength asymmetry" = mean(as.matrix(depHprime - 
                                                                depLprime), na.rm = TRUE)
      }
    }
    if ("SA" %in% index) {
      di <- dfun(web)$dprime
      dj <- dfun(t(web))$dprime
      if (SAmethod == "log") {
        lgmeani <- mean(log(di[di > 0]))
        lgmeanj <- mean(log(dj[dj > 0]))
        SA <- (lgmeanj - lgmeani)/sum(lgmeani, lgmeanj)
      }
      if (SAmethod == "Bluethgen") {
        wmeani <- sum(di * rowSums(web.e))/sum(web.e)
        wmeanj <- sum(dj * colSums(web.e))/sum(web.e)
        SA <- (wmeanj - wmeani)/sum(wmeani, wmeanj)
      }
      out$"specialisation asymmetry" <- SA
    }
    if (any(c("linkage density", "weighted connectance") %in% 
            index)) {
      preytot.mat <- matrix(rep(colSums(web), NROW(web)), 
                            NROW(web), byrow = TRUE)
      preyprop.mat <- web/preytot.mat
      predtot.mat <- matrix(rep(rowSums(web), NCOL(web)), 
                            NROW(web), byrow = FALSE)
      predprop.mat <- web/predtot.mat
      if (logbase == 2 | logbase == "2") {
        H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x * 
                                                          log2(x), na.rm = TRUE))
        H_Pk <- apply(predprop.mat, 1, function(x) -sum(x * 
                                                          log2(x), na.rm = TRUE))
        n_Nk <- ifelse(colSums(web) != 0, 2^H_Nk, 0)
        n_Pk <- ifelse(rowSums(web) != 0, 2^H_Pk, 0)
      }
      if (logbase == "e") {
        H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x * 
                                                          log(x), na.rm = TRUE))
        H_Pk <- apply(predprop.mat, 1, function(x) -sum(x * 
                                                          log(x), na.rm = TRUE))
        n_Nk <- ifelse(colSums(web) != 0, exp(H_Nk), 
                       0)
        n_Pk <- ifelse(rowSums(web) != 0, exp(H_Pk), 
                       0)
      }
      V <- sum(rowSums(web)/sum(web) * n_Pk)
      G <- sum(colSums(web)/sum(web) * n_Nk)
      LD_q <- 0.5 * (V + G)
      if ("linkage density" %in% index) 
        out$"linkage density" <- LD_q
      if ("weighted connectance" %in% index) 
        out$"weighted connectance" <- LD_q/sum(dim(web))
    }
    if ("Fisher alpha" %in% index) {
      fish <- try(fisherfit(web)$estimate, silent = TRUE)
      if (inherits(fish, "try-error")) {
        out$"Fisher alpha" <- NA
      }
      else {
        out$"Fisher alpha" <- fish
      }
    }
    if (any(c("interaction evenness", "Alatalo interaction evenness", 
              "Shannon diversity") %in% index)) {
      p_i.mat <- web/sum(web)
      SH <- -sum(p_i.mat * log(p_i.mat), na.rm = TRUE)
      if ("Shannon diversity" %in% index) 
        out$"Shannon diversity" <- SH
      IE <- ifelse(intereven == "prod", SH/log(prod(dim(web))), 
                   SH/log(sum(web > 0)))
      if ("interaction evenness" %in% index) 
        out$"interaction evenness" <- IE
      if ("Alatalo interaction evenness" %in% index) {
        evenness <- function(web) {
          pk <- web/sum(web)
          (Alatalo <- (1/sum(pk^2) - 1)/(exp(-sum(pk * 
                                                    log(pk), na.rm = TRUE)) - 1))
        }
        E <- evenness(web)
        out$"Alatalo interaction evenness" <- E
      }
    }
    if ("H2" %in% index) {
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - 
                                                                         round(x)) < tol
      if (any(is.wholenumber(web) == FALSE)) 
        H2_integer <- FALSE
      H2 <- as.numeric(H2fun(web, H2_integer = H2_integer)[1])
      out$H2 = ifelse(H2 < 0, 0, H2)
    }
    netw.index <- match(c("connectance", "web asymmetry", 
                          "links per species", "number of compartments", 
                          "compartment diversity", "nestedness", 
                          "NODF", "weighted nestedness", "weighted NODF", 
                          "ISA", "SA", "interaction evenness", 
                          "Alatalo interaction evenness", "Fisher alpha", 
                          "H2", "Shannon diversity", "linkage density", 
                          "weighted connectance"), index)
    exclude.index <- netw.index[!is.na(netw.index)]
    gindex <- if (length(exclude.index) == 0) 
      index
    else index[-exclude.index]
    if (length(gindex) > 0) 
      outg <- grouplevel(web, index = gindex, level = level, 
                         weighted = weighted, extinctmethod = extinctmethod, 
                         nrep = nrep, CCfun = CCfun, dist = dist, normalise = normalise, 
                         empty.web = empty.web, logbase = logbase, fcweighted = fcweighted, 
                         fcdist = fcdist)
    if (exists("outg")) {
      if (is.list(outg)) {
        SEQ <- seq(1, 2 * length(outg[[1]]), by = 2)
        sorted.outg <- c(outg[[1]], outg[[2]])
        outg <- sorted.outg[order(c(SEQ, SEQ + 1))]
      }
      out <- c(unlist(out), outg)
    }
    else {
      out <- unlist(out)
    }
  }
  if (legacy == TRUE) {
    out <- .networklevel(web, index = index, ISAmethod = ISAmethod, 
                         SAmethod = SAmethod, extinctmethod = extinctmethod, 
                         nrep = nrep, plot.it.extinction = FALSE, plot.it.dd = FALSE, 
                         CCfun = CCfun, dist = dist, normalise = normalise, 
                         empty.web = empty.web, logbase = logbase, fcweighted = fcweighted, 
                         fcdist = fcdist)
  }
  return(out)
}

function (web, index = "ALLBUTDD", level = "both", 
          weighted = TRUE, empty.web = TRUE, dist = "horn", CCfun = mean, 
          logbase = "e", normalise = TRUE, extinctmethod = "r", 
          nrep = 100, fcdist = "euclidean", fcweighted = TRUE) 
{
  if (level == "both") 
    for.higher <- for.lower <- TRUE
  if (level == "higher") {
    for.higher <- TRUE
    for.lower = FALSE
  }
  if (level == "lower") {
    for.lower <- TRUE
    for.higher = FALSE
  }
  if ("vulnerability" %in% index & level == "higher") 
    warning("You requested 'vulnerability' for the higher level, although it is not a higher level index! You will get 'generality' instead (same thing, really).", 
            call. = FALSE)
  if ("generality" %in% index & level == "lower") 
    warning("You requested 'generality' for the lower level, although it is not a lower level index! You will get 'vulnerability' instead (same thing, really).", 
            call. = FALSE)
  if (for.higher) 
    outh <- one.grouplevel(web, level = "higher", index = index, 
                           weighted = weighted, empty.web = empty.web, dist = dist, 
                           CCfun = CCfun, logbase = logbase, normalise = normalise, 
                           extinctmethod = extinctmethod, nrep = nrep, fcdist = fcdist, 
                           fcweighted = fcweighted)
  if (for.lower) 
    outl <- one.grouplevel(web, level = "lower", index = index, 
                           weighted = weighted, empty.web = empty.web, dist = dist, 
                           CCfun = CCfun, logbase = logbase, normalise = normalise, 
                           extinctmethod = extinctmethod, nrep = nrep, fcdist = fcdist, 
                           fcweighted = fcweighted)
  if (level == "higher") 
    return(outh)
  if (level == "lower") 
    return(outl)
  if (!("degree distribution" %in% index)) {
    out <- c(outh, outl)
    SEQ <- seq(1, length(out), by = 2)
    out <- out[order(c(SEQ, SEQ + 1))]
  }
  if (level == "both") {
    if (!("degree distribution" %in% index)) {
      out <- c(outh, outl)
      SEQ <- seq(1, length(out), by = 2)
      out <- out[order(c(SEQ, SEQ + 1))]
    }
    else {
      out <- list(HL = outh, LL = outl)
    }
    return(out)
  }
}