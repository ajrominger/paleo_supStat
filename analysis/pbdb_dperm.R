# **script to caculate d stats on sstat objects and null permutations**

library(socorro)
R.utils::sourceDirectory('R', modifiedOnly = FALSE)

# read in data
pbdbGenDiv <- read.csv('data/pbdb_3TPub_genera.csv', as.is = TRUE)
pbdbTax <- read.csv('data/pbdb_taxa.csv', as.is = TRUE)
# pbdbFamDiv <- read.csv('data/pbdb_3TPub-corrected.csv', row.names = 1)
load('data/pbdb_sstat_objects.RData')

# the indeces of otu names in `pbdbTax` ordered by their occurence in `pbdbGenDiv`
# needed to match permuted families to genera
genHash <- match(pbdbGenDiv$otu, pbdbTax$otu)


#' function to calculate KS test d-stat on sstat objects
#' @param x an sstat object

ks.sstat <- function(x, ...) {
    # browser()
    dat <- abs(unlist(x$Px.sub))
    # dat <- dat
    # `PPx(x)` was designed to be plotted as inverse CDF on absolute values 
    # as in `plot.sstat`, so undo that with `1 - 0.5 * PPx`
    # pfun <- function(X) 1 - 0.5 * x$PPx(X)
    # pfun <- function(X) x$PPx(exp(X))
    # pfun <- function(X) 0.5 + 0.5 * x$PPx(X, comp = FALSE)
    pfun <- function(X, ...) log(x$PPx(X, ...))
    
    n <- length(dat)
    pobs <- log(rev((1:n) / n))
    out <- pfun(sort(dat), comp = TRUE, ...) #- pobs
    
    plot(sort(dat), pobs)
    points(sort(dat), out, col = 'red')
    range(out - pobs)

     
    # plot(simpECDF(dat), type = 'l')
    # lines(sort(dat), pfun(sort(dat)), col = 'red')
    
    # plot(out)
    # points(1 / n - out, col = 'red')
    
    return(max(out, 1 / n - out))
}


foo <- sstatPBDBPhy
ks.sstat(foo)


plot(simpECDF(unlist(foo$Px.sub)))


# sstat on real (non-permuted) data
dObsFam <- ks.sstat(sstatPBDBfam3TP)
dObsOrd <- ks.sstat(sstatPBDBOrd)
dObsCls <- ks.sstat(sstatPBDBCls)
dObsPhy <- ks.sstat(sstatPBDBPhy)



n <- length(unlist(sstatPBDBPhy$Px.sub))
foo <- max(abs((1 - 0.5 * sstatPBDBfam3TP$PPx(sort(unlist(sstatPBDBfam3TP$Px.sub))) - (1:n) / n)))


plot(1 - 0.5 * sstatPBDBPhy$PPx(sort(unlist(sstatPBDBPhy$Px.sub))), (1:n) / n)
abline(0, 1, col = 'red')

plot(simpECDF(unlist(sstatPBDBfam3TP$Px.sub)))
curve(1 - 0.5 * sstatPBDBfam3TP$PPx(x), col = 'red', add = TRUE)
segments(0.7, 0.7, 0.7, 0.7-foo)

curve(0.5 + 0.5 * sstatPBDBfam3TP$PPx(x), from = -10, to = 10)


# repeatedly permute data and calculate null ks statistics
B <- 999
dNull <- parallel::mclapply(1:B, mc.cores = 3, FUN = function(i) {
    newFam <- sample(pbdbTax$family)
    newDiv <- tidy2mat(pbdbGenDiv$tbin, newFam[genHash], pbdbGenDiv$T3PubDiv)
    newFlux <- calcFlux(newDiv)
    newSstat <- sstatComp(newFlux, minN = 10, plotit = FALSE)
    
    return(ks.sstat(newSstat))
})

dNull <- c(unlist(dNull), dObs)


save(dNull, file = 'dnull_temp.RData')
