## 6: Plotting

## Combines the node-level and M-H distances into a single dataframe, then plots
## the node-level and network-level violin plots


pb_download("zscore50_5.RData",
            dest="data",
            tag="data.v.1")
pb_download("diffDist5Zscore.RData",
            dest="data",
            tag="data.v.1")
pb_download("sexlvlProp50Z.RData",
            dest="data",
            tag="data.v.1")
load("data/zscore50_5.RData")
load("data/diffDist5Zscore.RData")
load("data/sexlvlProp50Z.RData")
load("data/nets_mix_clean.RData")
load('data/spec_all.RData')

library(reshape2)
library(ggplot2)
library(stringr)
library(tidyverse)

##----------------
#Node-level plot:

##combine network statistics and M-H distance in a single dataframe. 
#Z-scores were used here to better show the magnitude of differences graphically
diff_traits5Zscore.df<-traitFrame(zscore50_5.df,
                                  diffDist5Zscore)

#Melt data frame into a form that ggplot likes better
nodelvlmeltZ<-melt(diff_traits5Zscore.df[,c(1,4,2,3,6,5)],
                   ID="SpSiteYr")
  #This orders the factors in the plot to match the order they're presented in the text


#Node-level plot
ggplot(nodelvlmeltZ,aes(x=variable,y=value))+
  geom_violin(draw_quantiles=0.5) +
  #geom_jitter(size=1,width=.2,shape=21) + #Alternative to add actual datapoints to the graph
  geom_hline(yintercept=c(0,1.96,-1.96),linetype=c("solid","dashed","dashed"),color=c("red","red","red"))+
  theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14))+
  ylab("Z-score Value")+
  xlab("Network Metric")+
  scale_x_discrete(labels=str_wrap(c("Degree", 
                                     "d'",
                                     "Weighted betweeness",
                                     "Weighted closeness",
                                     "M-H distance"),
                                   width = 10))




nodelvlmeltZ1<-melt(diff_traits5Zscore.df[,c(1,5,3,4,10,6,8)],
                   ID="SpSiteYr")

ggplot(nodelvlmeltZ1,aes(x=variable,y=value))+
  geom_violin(draw_quantiles=0.5) +
  #geom_jitter(size=1,width=.2,shape=21) + #Alternative to add actual datapoints to the graph
  geom_hline(yintercept=c(0,1.96,-1.96),linetype=c("solid","dashed","dashed"),color=c("red","red","red"))+
  theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14))+
  ylab("Z-score Value")+
  xlab("Network Metric")+
  scale_x_discrete(labels=str_wrap(c("Degree", 
                                     "d'", 
                                     "Weighted betweeness",
                                     "Weighted closeness",
                                     "M-H distance"),
                                   width = 10))


##-------------
## Network-level plot

netMetZ <- cbind(sexlvlProp50Z,robZ[,c(1,2)])
sexlvlmeltZ<-melt(sexlvlProp50Z[,c(3,4,2,1,5,6,7,8)],
                  ID="SpSiteYr")
  #This orders the factors in the plot to match the order they're presented in the text

ggplot(sexlvlmeltZ,aes(x=variable,y=value))+
  geom_violin(draw_quantiles=0.5)+
  #geom_jitter(size=1,width=.2,shape=21,color= "gray65") + # to add the actual points to the plot
  geom_hline(yintercept=c(0,1.96,-1.96),linetype=c("solid","dashed","dashed"),color=c("red3","red3","red3"),size =0.9)+
  ylim(-12.5,12.5)+
  ylab("Z-score Value")+
  xlab("Network Metric")+
  theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14))+
  scale_x_discrete(labels=str_wrap(c("Pollinator niche overlap",
                                     "Plant niche overlap",
                                     "H2'",
                                     "NODF", 
                                     "Pollinator redundancy",
                                     "Plant redundancy",
                                     "Pollinator robustness",
                                     "Plant robustness"),
                                   width = 13))



##-----Summary of dataset--------

#total networks in dataset
load('data/spec_all.RData')
length(unique(spec_all$SiteYr))
  #312


#total networks used in analyses (i.e., had at least 1 species where both sexes present)
load('data/nets_mix_clean2k.RData')
length(nets_mix_clean2k[[1]])
  #256


#total number of pollinator species
net_list <- sub('\\.',' ',names(nets_mix_clean[[1]]))
spec_all_analysis <- spec_all[spec_all$SiteYr %in% net_list,]
length(unique(spec_all_analysis$GenusSpecies))
  #393

#total number of plant species
length(unique(spec_all_analysis$PlantGenusSpecies))
  #260

#unique interactions
length(unique(paste(spec_all_analysis$PlantGenusSpecies,
                    spec_all_analysis$GenusSpecies)
              ))
  #2968

##--Supplementary table: Species presence in different datasets

#list of sites with at least 1 species having both sexes
finalsites.ls <- sub("\\.", " ",names(nets_mix_clean2k[[1]]))

spec_final <- spec_all[spec_all$SiteYr %in% finalsites.ls,]

networkSums <- function(data, data_location) {
  spec_data <- filter(spec_final,dataset==data_location)
  
  loc.ls <- lapply(unique(spec_data$GenusSpecies), function(x){
    sp <- spec_data[spec_data$GenusSpecies==x,]
    num <- length(unique(sp$SiteYr))
    df <- as.data.frame(x)
    df$y <- num
    names(df) <- c("Species",data_location)
    return(df)
  })
  
  loc.df <- do.call(rbind,loc.ls)
  loc.df$Family <- spec_final$Family[match(loc.df$Species,spec_final$GenusSpecies)]
  return(loc.df)
}

#number of networks containing each species in each dataset
yos.df <- networkSums(spec_final,"y")
hr.df <- networkSums(spec_final,"h")
si.df <- networkSums(spec_final,"s")

#merging datasets
all.df <- full_join(hr.df,yos.df,by="Species")
all.df <- full_join(all.df,si.df,by="Species")

##adding which species were retained in the node-level analyses
sp.df <- as.data.frame(trimws(unique(sexDiffs5.df[[1]]$species)))
sp.df$Retained <- rep(1)
names(sp.df) <- c("Species","Retained")
all.df <- left_join(all.df,sp.df,by="Species")


#cleaning up, 
all.df <- all.df[,c(3,1,2,4,6,8)]
names(all.df) <- c("Family","Species","HR","YOS","SI","Retained")
all.df$Family <- spec_all$Family[match(all.df$Species,spec_all$GenusSpecies)]

write.csv(all.df,"data/NetworkNumbers.csv")


