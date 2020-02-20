## 4: Testing and models


pb_download("zscore50_2.RData",
            dest="data",
            tag="data.v.1")

pb_download("sexDiffsProp50YH_20.Rdata",
            dest="data",
            tag="data.v.1")

load("data/zscore50_2.RData")
load("data/sexDiffsProp50YH_20.Rdata")

metric.ls <- c("degree","species.strength","weighted.betweenness",
               "weighted.closeness","d")


###------------------
##Test: proportion of species+sites where m v f difference in
##observed network was larger than many of the simulations

overallTest(sexDiffsProp50_2.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_3.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_5.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_10.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_20.df, metric.ls, tails=1, zscore=F)



overallTest(zscore50_2.df, metric.ls, zscore=T)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=2)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=-1)


#splevel
spLevelTest(zscore50_2.df,metric.ls,zscore=T)
#results: tons of zeros

#familylevel
spLevelTest(sexDiffsProp50_5.df,metric.ls,zscore=F,level="Family")


sexDiffsProp50_2.df %>%
  mutate(Sp = gsub( "_.*$", "", SpSiteYr))%>%
  arrange(desc(Sp)) -> sexDiffsProp50_2.df

zscore50_2.df %>%
  mutate(Sp = gsub( "_.*$", "", SpSiteYr))





##density plot for degree, absolute differences
plot(density(zscore50_2.df$degree, na.rm = T))
abline(v=0)
abline(v=1.645, lty="dashed")
abline(v=-1.645, lty="dashed")

plot(density(zscore50_2.df$species.strength, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$weighted.betweenness, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$weighted.closeness, na.rm = T))
abline(v=0)

plot(density(zscore50_2.df$d, na.rm = T))
abline(v=0)


#regress amt of difference against absolute specialization?
#i.e, are more generalized sp more different b/w males and females?






##
#----------Roswell paper notes

#6 meadow sites, sampling at minimum 6 30-mintue sampling periods
  #for 3 days
  #5 rounds throughout the summer each

#Do males and females overlap in diet?
  #permuted bee sex in the visitation record
  #9999 iterations: this many reqed to stabilize for 0.05 alpha
    #based on North, Curtis, and Sham 2002
  #)"When the observed dissimilarity was greater than 9500
    #of the 9999 simulated dissimilarities, we concluded that
    #we had detected a difference in the pattern of floral
    #visitation between conspecific male and female bees, 
    #given the observed diet breadth and abundance of each sex.
  #dissimilarity calculated using: Morisita-Horn index
    #good for mixed-size and small sample data
  #comparison at the SPECIES level, across all sites and 
    #sample rounds (several w/in 2016)




