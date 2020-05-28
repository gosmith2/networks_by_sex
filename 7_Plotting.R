## 7: Plotting

## Combines the node-level and M-H distances into a single dataframe, then plots
## the node-level and network-level violin plots

load("data/zscore50_2.RData")
load("data/diffDistZscore.RData")
load("daata/sexlevelProp50Z.RData")


##----------------
#Node-level plot:

##combine network statistics and M-H distance in a single dataframe. 
#Z-scores were used here to better show the magnitude of differences graphically
diff_traits5Zscore.df<-traitFrame(zscore50_5.df,diffDist5Zscore)


#Melt data frame into a form that ggplot likes better
nodelvlmeltZ<-melt(diff_traits5Zscore.df[,c(1,2,3,4,5,6,8)],ID="SpSiteYr")


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

sexlvlmeltZ<-melt(sexlvlProp50Z[,4:14],ID="SpSiteYr")
sexlvlmeltzSm <- filter(sexlvlmeltZ,variable!="functional.complementarity.HL",variable!="functional.complementarity.LL")
ggplot(sexlvlmeltzSm,aes(x=variable,y=value))+
  geom_violin(draw_quantiles=0.5)+
  # geom_jitter(size=1,width=.2,shape=21,color= "gray65") +
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




