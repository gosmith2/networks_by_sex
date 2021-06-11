## 6: Plotting

## Combines the node-level and M-H distances into a single dataframe, then plots
## the node-level and network-level violin plots


pb_download("zscore50_O.RData",
            dest="data",
            tag="data.v.1")
pb_download("diffDist5ZscoreO.RData",
            dest="data",
            tag="data.v.1")
pb_download("sexlvlProp50Z.RData",
            dest="data",
            tag="data.v.1")
load("data/zscore50_5.RData")
load("data/diffDist5Zscoreboot.RData")
load("data/sexlvlProp50Z.RData")
load("data/nets_mix_clean.RData")
load('data/spec_all.RData')

library(reshape2)
library(ggplot2)
library(stringr)
library(tidyverse)

##----------------
#Node-level plot:


#Melt data frame into a form that ggplot likes better, re-order columns to match text
nodelvlmeltO<-melt(zjoinO[,c(1,3,2,5,4)],
                   ID="SpSiteYr")
nodelvlmeltS<-melt(zjoinS[,c(1,3,2,5,4)],
                   ID="SpSiteYr")


#Node-level plot
ggplot(nodelvlmeltS,aes(x=variable,y=value))+
  geom_violin(draw_quantiles=0.5,size=0.9,aes(fill=variable),show.legend = F) +
  scale_fill_viridis(discrete=T,alpha=0.5)+
  #geom_jitter(size=1,width=.2,shape=21) + #Alternative to add actual datapoints to the graph
  geom_hline(yintercept=0,linetype="dashed",color="red",size=0.9)+
  theme(panel.background = NULL,
        axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_blank())+
  ylab("Z-score Value")+
  xlab("Network Metric")+
  scale_x_discrete(labels=str_wrap(c("Degree", 
                                     "d'",
                                     "Weighted closeness",
                                     "M-H distance"),
                                   width = 10))



##-------------
## Network-level plot

netMeltO <- melt(zNetO[,c(2,1,3,4,5,6,7)],ID='SpSiteYr')
netMeltS <- melt(zNetS[,c(2,1,3,4,5,6,7)],ID='SpSiteYr')


ggplot(netMeltS,aes(x=variable,y=value,fill=variable,alpha=0.5))+
    theme(panel.background = NULL)+
    geom_violin(draw_quantiles=0.5,size=0.9,show.legend = F) +
    scale_fill_manual(values=c("#51127CFF","#B63679FF","#FB8861FF","#FB8861FF","#FCFDBFFF","#FCFDBFFF"))+
    #scale_fill_viridis(discrete=T,alpha=0.5,option="A")+
    #geom_jitter(size=1,width=.2,shape=21,color= "gray65") + # to add the actual points to the plot
    geom_hline(yintercept=0,linetype="dashed",color="red3",size =0.9)+
    ylim(-12.5,12.5)+
    ylab("Z-score Value")+
    xlab("Network Metric")+
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          axis.text.y=element_text(size=12),
          axis.title=element_text(size=14))+
  scale_x_discrete(labels=str_wrap(c("H2'",
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
length(nets_mix_clean10kO[[1]])
  #256


#total number of pollinator species
net_list <- sub('\\.',' ',names(nets_mix_clean10kO[[1]]))
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
finalsites.ls <- sub("\\.", " ",names(nets_mix_clean10k[[1]]))

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

