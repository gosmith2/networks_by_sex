

samp2site.spp <- function(site, spp, abund, FUN=mean) {
  ## This functions takes site-species-abundance data and creates a
  ## matrix where the sites are columns and the rows are species.
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}



fix.white.space <- function(d) {
  ## function to clean up white-space in a column of data (replaces all
  ## instances of white-space with " " and empty cells with ""
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
