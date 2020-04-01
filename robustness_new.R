
simExtinction <- function(nets,
                          extinction.method,
                          spec,
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
