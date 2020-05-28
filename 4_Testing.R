## 4: Overall test of network role differences

#Generates a test to summarize the number of networks where sex
#associated differences are larger than the null expectation. 

#Setup: Download and attach data, specify metrics to test
pb_download("zscore50_5.RData",
            dest="data",
            tag="data.v.1")

pb_download("sexDiffsProp50_5.Rdata",
            dest="data",
            tag="data.v.1")

load("data/zscore50_5.RData")
load("data/sexDiffsProp50_5.Rdata")

metric.ls <- c("degree","species.strength","weighted.betweenness",
               "weighted.closeness","d")


###------------------
##Test: proportion of species+sites where m v f difference in
##observed network was larger than many of the simulations. The output
##is the proportion of observations (Sp+Site+Year) where the male-female
## difference diverged from null expectations more than some threshold 
## (with the specific threshold used based on the tails of the test)


overallTest(sexDiffsProp50_5.df, metric.ls, tails=2, zscore=F)





######-------------------------------

# - Alternative code that will be removed before pub

######-------------------------------

overallTest(sexDiffsProp50_2.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_3.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_10.df, metric.ls, tails=1, zscore=F)
overallTest(sexDiffsProp50_20.df, metric.ls, tails=1, zscore=F)




overallTest(zscore50_2.df, metric.ls, zscore=T)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=2)

overallTest(zscore50_2.df, metric.ls, zscore=T, tails=-1)


#splevel
spLevelTest(zscore50_5.df,metric.ls,zscore=T)
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




