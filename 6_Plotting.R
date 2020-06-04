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


##----------------
#Node-level plot:

##combine network statistics and M-H distance in a single dataframe. 
#Z-scores were used here to better show the magnitude of differences graphically
diff_traits5Zscore.df<-traitFrame(zscore50_5.df,
                                  diffDist5Zscore)


#Melt data frame into a form that ggplot likes better
nodelvlmeltZ<-melt(diff_traits5Zscore.df[,c(1,2,3,4,5,6,8)],
                   ID="SpSiteYr")


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
                                     "Species strength", 
                                     "Weighted betweeness",
                                     "Weighted closeness",
                                     "d'",
                                     "M-H distance"),
                                   width = 10))


##-------------
## Network-level plot

sexlvlmeltZ<-melt(sexlvlProp50Z[,1:8],
                  ID="SpSiteYr")

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
  scale_x_discrete(labels=str_wrap(c("NODF", 
                                     "H2'", 
                                     "Pollinator niche overlap",
                                     "Plant niche overlap",
                                     "Pollinator robustness",
                                     "Plant robustness",
                                     "Generality",
                                     "Vulnerability"),
                                   width = 13))



##-----Summary of dataset--------

#total networks in dataset
load('data/spec_all.RData')
length(unique(spec_all$SiteYr))
  #312


#total networks used in analyses (i.e., had at least 1 species where both sexes present)
load('data/nets_mix_clean.RData')
length(nets_mix_clean[[1]])
  #257


#total number of pollinator species
net_list <- sub('\\.',' ',names(nets_mix_clean[[1]]))
spec_all_analysis <- spec_all[spec_all$SiteYr %in% net_list,]
length(unique(spec_all_analysis$GenusSpecies))
  #422

#total number of plant species
length(unique(spec_all_analysis$PlantGenusSpecies))
  #267

#unique interactions
length(unique(paste(spec_all_analysis$PlantGenusSpecies,
                    spec_all_analysis$GenusSpecies)
              ))
  #3044
