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

breakNetSex <- function(spec.dat, site, year, mix){
  ## puts data together in a list and removes empty matrices
  #browser()
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

#below is duplicate!!!!!!!!!!!!!
breakNetSex <- function(spec.dat, site, year, mix){
  ## puts data together in a list and removes empty matrices
  agg.spec <- aggregate(list(abund=spec.dat$GenusSpeciesSex),
                        list(GenusSpeciesSex=spec.dat$GenusSpeciesSex,
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
                  spp=y[,"GenusSpeciesSex"],
                  abund=y[,"abund"])
  })
  return(comms)
}
#above is duplicate


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

calcSpec <- function(nets, spec){
  ## applies specieslevel from bipartite to networks
  species.lev <- lapply(nets, function(x){
    sl <- specieslevel(x)
    return(sl)
  })
  
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

ran.sex<-function(spec.data){
  col.ls<-lapply(unique(spec.data$SiteYr),function(x){
    net<-filter(spec.data,SiteYr==x)
    sp.mix.col <-lapply(unique(net$GenusSpecies), function(y){
   #   browser()
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

#mclapply in the parallel package

ran.gen<-function(spec.data,iterations){
  #setup: add column to actual observation df and initiate list
  spec.data$MixSex<-spec.data$Sex
  spec.data$GenusSpeciesMix<-paste(spec.data$GenusSpecies,
                                   spec.data$MixSex,
                                   sep="_")
  destList<-list() 
  destList[[1]]<-spec.data
  
  #scramble sex column for within each species within each site
  it.vec<-1:iterations
  randoms<-mclapply(it.vec,function(z){
    spec.data<-ran.sex(spec.data)
    spec.data$GenusSpeciesMix<-paste(spec.data$GenusSpecies,
                                     spec.data$MixSex,
                                     sep="_")
    return(spec.data)
  },mc.core=6)
  destList<-c(destList,randoms)
  return(destList)
}
