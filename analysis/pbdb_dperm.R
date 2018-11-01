# **script to caculate d stats on sstat objects and null permutations**

library(socorro)
R.utils::sourceDirectory('R', modifiedOnly = FALSE)

# read in data
pbdbGenDiv <- read.csv('data/pbdb_3TPub_genera.csv', as.is = TRUE)
pbdbTax <- read.csv('data/pbdb_taxa.csv', as.is = TRUE)

# the indeces of otu names in `pbdbTax` ordered by their occurence in `pbdbGenDiv`
# needed to match permuted families to genera
genHash <- match(pbdbGenDiv$otu, pbdbTax$otu)


#' function to calculate KS test d-stat on sstat objects
#' @param x an sstat object

ks.sstat <- function(x, ...) {
    dat <- unlist(x$Px.sub)
    # `PPx(x)` was designed to be plotted as inverse CDF on absolute values 
    # as in `plot.sstat`, so undo that with `1 - 0.5 * PPx`
    pfun <- function(dat) 1 - 0.5 * x$PPx(dat) 
    
    n <- length(dat)
    out <- pfun(sort(dat), ...) - (0:(n - 1)) / n
    
    return(max(out, 1 / n - out))
}

# sstat on real (non-permuted) data
pbdbFamFlux <- calcFlux(pbdbFamDiv)
sstatPBDBfam3TP <- sstatComp(pbdbFamFlux, minN = 10, plotit = FALSE)
dObs <- ks.sstat(sstatPBDBfam3TP)

# repeatedly permute data and calculate null ks statistics
B <- 999
dNull <- parallel::mclapply(1:B, mc.cores = 8, FUN = function(i) {
    newFam <- sample(pbdbTax$family)
    newDiv <- tidy2mat(pbdbGenDiv$tbin, newFam[genHash], pbdbGenDiv$T3PubDiv)
    newFlux <- calcFlux(newDiv)
    newSstat <- sstatComp(newFlux, minN = 10, plotit = FALSE)
    
    return(ks.sstat(newSstat))
})

dNull <- c(unlist(dNull), dObs)
