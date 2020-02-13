## 5: Dissimilarity indexing

library(vegan)

load("data/mix_netsYHS.RData")

obs <- nets.mix.clean[[1]]

obsz<-obs[["Zamora.2014"]]

mh(as.matrix(obsz))

mh(obsz[2:3,c("Halictus tripartitus_f","Halictus tripartitus_m")]))

vegdist(t(obsz),method="horn") #(make this a matrix, rather than a distance object)



test1<-distComp(nets.mix.clean,"horn")

spec.h %>%
  filter(Site=="Berm",Year==2012,GenusSpecies=="Toxomerus marginatus")%>%
  select(GenusSpecies,Sex)
