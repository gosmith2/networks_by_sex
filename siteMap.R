

library(tidyverse)    # our old friend
library(rgdal)        # main spatial package
library(maptools)     # tools for manipulating spatial data
library(ggmap)        # for downloading maps 
library(rgeos)        # more map tools
library(worldmap)
library(measurements)

load("data/specimens-hr.RData",verbose=TRUE) 
#For whatever reason, resulting df is called "dd"
spec.h<-dd

spec.y <-read.csv("data/specimens-yos.csv")

spec.h$lat1 <- sub(" W.*$", "", spec.h$GPS)
spec.h$lat1 <- substring(spec.h$lat1,2)
spec.h$long1 <- sub("^[^W]*", "", spec.h$GPS)
spec.h$long1 <- substring(spec.h$long1,2)
spec.h$long <- measurements::conv_unit(spec.h$long1, from = 'deg_dec_min', to = "dec_deg")
spec.h$lat <- measurements::conv_unit(spec.h$lat1, from = 'deg_dec_min', to = "dec_deg")
head(spec.h$lat)

#nope, wont work yet. First, gotta replace all spaces with '. then 
$

keeps.spat <- c("Site",
                "Lat",
                "Long")


