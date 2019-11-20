

library(tidyverse)    # our old friend
library(rgdal)        # main spatial package
library(maptools)     # tools for manipulating spatial data
library(ggmap)        # for downloading maps 
library(rgeos)        # more map tools
library(rworldmap)
library(measurements)
library(ggspatial)

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
spec.h$Lat <- as.numeric(spec.h$lat)
spec.h$Long <- as.numeric(spec.h$long*-1)
head(spec.h$Long)

keeps.spat <- c("Site",
                "Lat",
                "Long")

bind_rows(select(spec.y,keeps.spat),
          select(spec.h,keeps.spat))->
  spec.all.spat

points.df <- data.frame(Lat = unique(spec.all.spat$Lat),
                        Long = unique(spec.all.spat$Long)
)
points.df <- points.df[1:78,]
points.df$Long <- points.df$Long*-1

minlat<-min(points.df$Lat)
maxlat<-max(points.df$Lat)
minlong<-min(points.df$Long)
maxlong<-max(points.df$Long)


plot(points.df)

newmap <- getMap(resolution = "high")
plot(newmap, ylim = c(37, 39), xlim = c(-124,-117), asp = 1)
points(points.df$Long,points.df$Lat, col="red")
#cool, so it works now. 

#next step: get in pasted on a google map

points.g <- SpatialPoints(points.df, proj4string=CRS("+init=epsg:4326"))

bbox.sites <- bbox(points.g)
#had to swap the order of these; think lat and long got confused somewhere?
bbox.all <- matrix(c(bbox.sites[2,1],
                     bbox.sites[1,1],
                     bbox.sites[2,2],
                     bbox.sites[1,2]),
                   ncol=2)


landscape <- ggmap(get_map(location=bbox.all, source="google", 
                           maptype="hybrid", crop=FALSE))

landscape + 
  geom_point(data= data.frame(points.g),
             aes(x=points.g$Long, y=points.g$Lat), color="red",cex=0.5)

