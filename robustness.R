


getExtinctionOrder <- function(veg.visit.degree,## visit/veg/degree
                               spec,
                               abund.col,
                               interaction.col,
                               order.extinct,
                               veg=NULL){
  ## function for simulating the extinction order of plants.
  ## if by visit, then abundaunce is calculated from networks,
  ## otherwise it is calculated from veg data.
  if(veg.visit.degree == "abundance"){
    abund <- aggregate(list(abund=spec[, abund.col]),
                       list(Year=spec$Year,
                            Site=spec$Site,
                            Species=spec[, abund.col]),
                       length)
    
    
  } else if(veg.visit.degree == "veg"){
    abund <- aggregate(list(abund =veg$logFlowerNum),
                       list(Year=veg$Year,
                            Site=veg$Site,
                            Species=veg[, abund.col]),
                       mean)
    
  } else if(veg.visit.degree == "degree"){
    abund <- aggregate(list(abund=spec[, interaction.col]),
                       list(Year=spec$Year,
                            Site=spec$Site,
                            Species=spec[, abund.col]),
                       function(x)length(unique(x)))
  }
  abund <- aggregate(list(abund=abund$abund),
                     list(Site=abund$Site,
                          Species=abund$Species),
                     mean)
  
  abund.dats <- split(abund, abund$Site)
  
  ord <- lapply(abund.dats, function(x){
    out <-  x[order(x$abund),]
    return(out)
  })
  gen.sp.order <- lapply(ord, function(x) x$Species)
  if(order.extinct == "high to low"){
    gen.sp.order <- lapply(gen.sp.order, rev)
  }
  
  return(gen.sp.order)
}

SimSecondExtinction <- function (web, participant = "higher",
                                 method = "abun", nrep = 10,
                                 details = FALSE, ext.row = NULL,
                                 ext.col = NULL){
  ## function modified from bipartite package to replicate
  ## extinction sequence because sveral species have the same
  ## abundance, leading to different estimates for robustness
  ## depending on the seed
  if (participant == "both" & method == "external")
    stop("Sorry, that won't work. When you specify the sequence, you have to choose one of the two levels. 'both' won't work.")
  if (!is.null(ext.row) & length(ext.row) != NROW(web))
    stop("The length of the external row vector is different from the numbers of rows in the network!")
  if (!is.null(ext.col) & length(ext.col) != NCOL(web))
    stop("The length of the external col vector is different from the numbers of cols in the network!")
  if (participant == "higher" & method == "external" & is.null(ext.col))
    stop("You need to provide an external sequence of extinction for the higher trophic level!")
  if (participant == "lower" & method == "external" & is.null(ext.row))
    stop("You need to provide an external sequence of extinction for the lower trophic level!")
  one.second.extinct <- function(web = web, participant = participant,
                                 method = method, ext.row = ext.row,
                                 ext.col = ext.col) {
    dead <- matrix(nrow = 0, ncol = 3)
    colnames(dead) <- c("no", "ext.lower", "ext.higher")
    m2 <- web
    i <- 1
    repeat {
      n <- sexExtinction(m2, participant = participant, method = method,
                      ext.row = ext.row, ext.col = ext.col)
      dead <- rbind(dead, c(i, attributes(m2 <- empty(n,
                                                      count = TRUE))$empty))
      if (participant == "lower" & NROW(m2) < 2)
        break
      if (participant == "higher" & NCOL(m2) < 2)
        break
      if (participant == "both" & min(dim(m2)) < 2)
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
      stop("PANIC!!")
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
  ## if (pmatch(method, c("abundance", "random", "degree", "external")) %in%
  ##     c(1, 3, 4)) {
  ##     out <- one.second.extinct(web = web, participant = participant,
  ##                               method = method, ext.row = ext.row, ext.col = ext.col)
  ## }
  o <- replicate(nrep, one.second.extinct(web = web, participant = participant,
                                          method = method,
                                          ext.row = ext.row, ext.col = ext.col),
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
  
  class(out) <- "bipartite"
  attr(out, "exterminated") <-
    c("both", "lower", "higher")[pmatch(participant,  c("both", "lower", "higher"))]
  out
}




simExtinction <- function(nets,
                          extinction.method,
                          spec,
                          participant,
                          ext.row.col,
                          nreps=10){ ## simulation repeats
  ## calculates the robustness of a network using Memmot et al.'s method
  ## takes the adjacency matrix, whether to drop species by abundance or
  ## degree and whther to drop the "higer" or "lower" level of the
  ## network
  ## returns a data frame with the site, robustness score
  sites <- sapply(strsplit(names(nets),  "\\."), function(x) x[[1]])
  ext <-  vector(mode="list", length=length(nets))
  browser()
  for(i in 1:length(nets)){
    if(all(dim(nets[[i]]) > 4)){
      if(participant == "lower"){
        net.names <- rownames(nets[[i]])
      }else if(participant == "higher"){
        net.names <- colnames(nets[[i]])
      }
      
      print(names(nets)[i])
      this.ext.order <- ext.row.col[[sites[i]]]
      this.ext.order <- this.ext.order[this.ext.order %in%
                                         net.names]
      rank.order <- match(net.names,
                          this.ext.order)
      if(participant == "lower"){
        this.ext <- SimSecondExtinction(nets[[i]],
                                        participant=participant,
                                        method=extinction.method,
                                        ext.row=rank.order,
                                        nrep=nreps)
      }else if(participant == "higher"){
        this.ext <- SimSecondExtinction(nets[[i]],
                                        participant=participant,
                                        method=extinction.method,
                                        ext.col=rank.order,
                                        nrep=nreps)
      }
    }else{
      this.ext <- NA
    }
    ext[[i]] <- this.ext
  }
  nulls <- sapply(ext, function(x) all(is.na(x)))
  net.names <-  names(nets)[!nulls]
  
  ext <- ext[!nulls]
  
  rob <- sapply(ext, robustness)
  
  split.names <- strsplit(net.names, "\\.")
  
  
  sites <- sapply(split.names,
                  function(x) x[1])
  years <-  sapply(split.names,
                   function(x) x[2])
  if(length(split.names[[1]]) > 2){
    SR <-  sapply(split.names,
                  function(x) x[3])
  } else {
    SR  <-  NA
  }
  dats <- data.frame(Site= sites,
                     Year=years,
                     SampleRound=SR,
                     Robustness=rob)
  rownames(dats) <- NULL
  return(dats)
}


sexExtinction <- function (web, participant = "both", method = "random", 
                           ext.row = NULL, ext.col = NULL) 
  #modified from bipartite so that if one sex of a given pollinator species goes extinct, the other sex will as well
{
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
    web <- web[sample(1:nrow(web)), sample(1:ncol(web))]
    rseq <- order(rowSums(web))
    cseq <- order(colSums(web))
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) {
      web[, cseq[1]] <- 0 
      ### look for any other columns with the same name; delete them too
      sp <- sub("_.*$", "", names(web[1,])[cseq[1]])
      web[,grepl(sp,names(web[1,]))] <- 0
    }
    if (partis.match == 3) {
      if (min(rowSums(web)) < min(colSums(web))) {
        web[rseq[1], ] <- 0
      }
      else {
        if (min(rowSums(web)) > min(colSums(web))) {
          web[, cseq[1]] <- 0
        }
        else {
          if (sample(2, 1) == 1) {
            web[rseq[1], ] <- 0
          }
          else {
            web[, cseq[1]] <- 0
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
      web[, ex] <- 0
    }
  }
  if (meths.match == 4) {
    rseq <- ext.row
    cseq <- ext.col
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) {
      web[, cseq[1]] <- 0
      ### zero out any other species with the same name as the removed species
      sp <- sub("_.*$", "", names(web[1,])[cseq[1]])
      web[,grepl(sp,names(web[1,]))] <- 0
    }
  }
  return(web)
}

ext.order <- getExtinctionOrder("abundance",
                                filter(spec.all,SiteYr=="Anderson 2010" | SiteYr=="BC1 2006"),
                                abund.col="GenusSpecies",
                                interaction.col="PlantGenusSpecies",
                                order.extinct=
                                  "low to high")

simExtinction(nets.mix.clean[[1]][1:2], extinction.method = "abundance", spec = filter(spec.all,SiteYr=="Anderson 2010" | SiteYr=="BC1 2006"),participant = "higher",ext.row.col = ext.order)
