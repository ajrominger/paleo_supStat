makePlot <- TRUE

oldwd <- setwd('~/Dropbox/Research/paleo_supStat')

# convenience function to produce a matrix of time by ord with cells
# of corrected diversity
source('code/pbdb_3t_pub.R')

# load other needed funcitons
paleoPlot <- function(...) {
    plot(..., xaxt = 'n')
    socorro::paleoAxis(1)
}
# source('~/R_functions/paleoPlot.R')

samp2site.spp <- socorro::tidy2mat
# source('~/R_functions/samp2site_spp.R')

logPlot <- function(..., log) {
    plot(..., log = log, 
         xaxt = ifelse(grepl('x', log), 'n', 's'), 
         yaxt = ifelse(grepl('y', log), 'n', 's'))
    s <- (1:2)[c(grepl('x', log), grepl('y', log))]
    
    socorro::logAxis(s)
}
# source('~/R_functions/logPlot.R')

my.ecdf <- socorro::simpECDF
# source('~/R_functions/my_ecdf.R')

source('code/sstat_comp.R')
source('code/sstat_methods.R')
source('code/Px_gam.R')

##########  load data  ##########
# setwd('data/old/pbdb_2013-05-28')

# raw occurence data
# pbdb.dat <- read.csv('marInv-occs.csv', as.is = TRUE)

pbdbNew <- read.csv('data/pbdb_data.csv', as.is = TRUE)

## convert column names
# names(pbdbNew)[names(pbdbNew) == 'tbin'] <- 'collections.10_my_bin'
# names(pbdbNew)[names(pbdbNew) == 'family'] <- 'occurrences.order_name'
# names(pbdbNew)[names(pbdbNew) == 'genus'] <- 'occurrences.genus_name'
# names(pbdbNew)[names(pbdbNew) == 'reference_no'] <- 'collections.reference_no'

# make a column for the million yr midpoint
pbdbNew$ma_mid <- (pbdbNew$max_ma + pbdbNew$min_ma) / 2

## use new instead of old
pbdb.dat <- pbdbNew


## get rid of poor temporal resolution
pbdb.dat <- pbdb.dat[pbdb.dat$tbin != '', ]

# get rid of bad taxonomy
pbdb.dat <- pbdb.dat[pbdb.dat$family != '', ]
# pbdb.dat <- pbdb.dat[pbdb.dat$occurrences.order_name != 'Ammonitida', ]
pbdb.dat <- pbdb.dat[pbdb.dat$genus != '', ]

# get bin times
pbdb.time <- sort(tapply(pbdb.dat$ma_mid, pbdb.dat$tbin, mean))
pbdb.dat$tbin <- factor(pbdb.dat$tbin, 
                                         levels = names(pbdb.time))


# data.frame of publication, diversity and 3T stat
ord.tbin.bias <- aggregate(list(div=pbdb.dat$genus),
                           list(ord=pbdb.dat$family,
                                tbin=pbdb.dat$tbin),
                           function(x) length(unique(x)))


## three timer stat

## matrix to determine three timers and part timers (sensu alroy 2008)
mt <- matrix(0, nrow = nlevels(pbdb.dat$tbin), 
             ncol = nlevels(pbdb.dat$tbin))
diag(mt) <- -10
mt[abs(row(mt) - col(mt)) == 1] <- 1

## loop through and compute three timers and part timers
timers <- lapply(split(pbdb.dat$tbin, pbdb.dat$genus), 
                 function(x) {
                     # browser()
                     tbins <- integer(nlevels(x))
                     tbins[as.integer(unique(x))] <- 1
                     t3 <- as.integer(mt %*% tbins == 2)
                     tp <- as.integer(mt %*% tbins == -8)
                     
                     return(cbind(t3, tp))
                 })

timers <- array(unlist(timers), dim = c(nrow(timers[[1]]), 2, length(timers)))

t3stat <- 1 - rowSums(timers[, 1, ]) / (rowSums(timers[, 1, ]) + rowSums(timers[, 2, ]))

ord.tbin.bias$T3.stat <- t3stat[match(ord.tbin.bias$tbin, 
                                      levels(pbdb.dat$tbin))]
ord.tbin.bias$T3.div <- ord.tbin.bias$div/ord.tbin.bias$T3.stat

# record pubs per tbin
tbin.pub <- tapply(pbdb.dat$reference_no,pbdb.dat$tbin,
                   function(x) length(unique(x)))
ord.tbin.bias$tbin.pub <- tbin.pub[ord.tbin.bias$tbin]


# calculate corrected diversity

pbdb.ord.div <- with(ord.tbin.bias,
                     pbdb.3t.pub(div, T3.stat, tbin.pub, ord, tbin, pbdb.time, 
                                 min.pub = 10, plotit = makePlot))

# corrected flux
pbdb.ord.flux <- apply(pbdb.ord.div, 2, function(x) {
    raw.flux <- diff(c(0, x))
    return(raw.flux[raw.flux != 0])
})

# sstat analysis
pbdb.sstat.ord.cor <- sstat.comp(pbdb.ord.flux, minN = 10, plotit = makePlot)
