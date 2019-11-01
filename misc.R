## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species.

samp2site.spp <- function(site, spp, abund, FUN=mean) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}

## does the reverse of samp2site
comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(Site = temp[, 1], Samp = temp[, 3],
             Date = temp[, 2])
}


## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
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

build.nets <- function(spec.dat,title) {
  SY <- breakNet(spec.dat, 'Site', 'Year')
  SSY <- breakNetSex(spec.dat,'Site','Year')
  GSY <- lapply(SY, graph.incidence, weighted=TRUE)
  GSSY <- lapply(SSY, graph.incidence, weighted=TRUE) 
  assign(paste(title,"SY",sep="_"), SY, envir = globalenv())
  assign(paste(title,"SSY",sep="_"), SSY, envir = globalenv())
  assign(paste(title,"GSY",sep="_"), GSY, envir = globalenv())
  assign(paste(title,"GSSY",sep="_"), GSSY, envir = globalenv())
  }


## return sorted unique values
id <- function(x) unique(sort(x))

## function to make pollinator visitation matrices
make.mats <- function(pollinator.id,
                      null.mat,
                      pollinator,
                      var1,
                      var2,
                      occ) {
  make.mat <- function(P) {
    var1 <- var1[pollinator==P]
    var2 <- var2[pollinator==P]
    occ <- occ[pollinator==P]
    null.mat[unique(var1),
             unique(as.character(var2))][!is.na(null.mat)] <- occ
    null.mat
  }
  mats <- lapply(pollinator.id, function(x) make.mat(x))
  names(mats) <- pollinator.id
  mats
}