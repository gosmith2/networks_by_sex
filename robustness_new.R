
calcRobCI <- function (dat.web, N, cLevels = c(0.95,0.9), participant, extinction.method) {
    ## calculate confidence intervals
    ci <- function(stats, nnull) {
        mean <- mean(stats, na.rm = T)
        sd <- sd(stats, na.rm = T)
        intervals <- lapply(cLevels,function(x){
            conf <- qt(x,df=nnull)*sd/sqrt(nnull+1)
            low <- mean - conf
            high <- mean + conf
            df <- data.frame(low,high)
            colnames(df) <- c(paste(x, "low", sep = "_"),paste(x, "high", sep = "_"))
            return(df)
        })
        
        interval.df <- do.call(cbind,intervals)
        return(interval.df)
    }
    
    ## check that matrix is proper format (no empty row/col and no NAs)
    if(all(is.na(dat.web) == FALSE)) {
        ## drop empty rows and columns
        dat.web <- as.matrix(empty(dat.web))
        ## check to make sure emptied matrix is large enough
        ## to calculate statistics on
        if(is.matrix(dat.web)){
            if(all(dim(dat.web) >= 2)) {
                ## calculate null metrics
                null.stat <- replicate(N,
                                       calcNullRob(dat.web,
                                                    null.fun= vaznull.fast,
                                                   extinction.method= extinction.method,
                                                   participant = participant),
                                       simplify=TRUE)
                ## calculate metrics from data
                true.stat <- calcMetric(dat.web,
                                        index=index)
                out.mets <- cbind(true.stat, null.stat)
                rowSeq <- seq(1:nrow(out.mets))
                ciValues <- lapply(rowSeq, function(x){
                    row <- out.mets[x,]
                    rowCI <- ci(row,N)
                    return(rowCI)
                })
                ciValues <- do.call(rbind,ciValues)
                ciValues$true <- true.stat 
                rownames(ciValues) <- names(true.stat)
                return(ciValues)
            }
        }
    }
    return(matrix(nrow=length(index)+4,ncol=(2*length(cLevels))+1))
}

calcNullRob <- function(dat.web,
                         null.fun,extinction.method = extinction.method,participant=participant,...) {
    sim.web <- null.fun(dat.web)
    return(simExtinction(sim.web,extinction.method=extinction.method,participant=participant,...))
}


simExtinction <- function(nets,
                          extinction.method,
                          participant="lower"){
    ## calculates the robustness of a network using Memmot et al.'s method
    ## takes the adjacency matrix, whether to drop species by abundance or
    ## degree and whther to drop the "higer" or "lower" level of the
    ## network
    ## returns a data frame with the site, robustness score and YPR
    ext <- lapply(nets, second.extinct,
                  participant=participant,
                  method=extinction.method)

    rob <- sapply(ext, robustness)

    sites <- sapply(strsplit(names(rob), "[.]"), function(x) x[1])
    years <- sapply(strsplit(names(rob), "[.]"), function(x) x[2])

    dats <- data.frame(Site= sites,
                       Year=years,
                       Robustness=rob)
    rownames(dats) <- NULL

    return(dats)
}


dropRandomSex <- function(data) {
    
}


