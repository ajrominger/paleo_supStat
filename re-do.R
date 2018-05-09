#######  corrects diversity by 3-timer then number of pubs  #######
makePlot <- TRUE


oldwd <- setwd('~/Dropbox/Research/paleo_supStat')

##  convenience function to produce a matrix of time by ord with cells
##  of corrected diversity
source('code/pbdb_3t_pub.R')

##	load other needed funcitons
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

##########  load data  ##########
setwd('data/old/pbdb_2013-05-28')

##	raw occurence data
pbdb.dat <- read.csv('marInv-occs.csv', as.is = TRUE)

## get rid of poor temporal resolution
pbdb.dat <- pbdb.dat[pbdb.dat$collections.10_my_bin != '', ]

##	get rid of bad taxonomy
pbdb.dat <- pbdb.dat[pbdb.dat$occurrences.order_name != '', ]
pbdb.dat <- pbdb.dat[pbdb.dat$occurrences.order_name != 'Ammonitida', ]
pbdb.dat <- pbdb.dat[pbdb.dat$occurrences.genus_name != '', ]

##	get bin times
pbdb.time <- sort(tapply(pbdb.dat$ma_mid, pbdb.dat$collections.10_my_bin, mean))
pbdb.dat$collections.10_my_bin <- factor(pbdb.dat$collections.10_my_bin, 
                                         levels = names(pbdb.time))


##  data.frame of publication, diversity and 3T stat
ord.tbin.bias <- aggregate(list(div=pbdb.dat$occurrences.genus_name),
                           list(ord=pbdb.dat$occurrences.order_name,
                                tbin=pbdb.dat$collections.10_my_bin),
                           function(x) length(unique(x)))


## three timer stat

## matrix to determine three timers and part timers (sensu alroy 2008)
mt <- matrix(0, nrow = nlevels(pbdb.dat$collections.10_my_bin), ncol = nlevels(pbdb.dat$collections.10_my_bin))
diag(mt) <- -10
mt[abs(row(mt) - col(mt)) == 1] <- 1

## loop through and compute three timers and part timers
timers <- lapply(split(pbdb.dat$collections.10_my_bin, pbdb.dat$occurrences.genus_name), 
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

ord.tbin.bias$T3.stat <- t3stat[match(ord.tbin.bias$tbin, levels(pbdb.dat$collections.10_my_bin))]
ord.tbin.bias$T3.div <- ord.tbin.bias$div/ord.tbin.bias$T3.stat

##	record pubs per tbin
tbin.pub <- tapply(pbdb.dat$collections.reference_no,pbdb.dat$collections.10_my_bin,function(x) length(unique(x)))
ord.tbin.bias$tbin.pub <- tbin.pub[ord.tbin.bias$tbin]


##	calculate corrected diversity

pbdb.ord.div <- with(ord.tbin.bias,
                     pbdb.3t.pub(div, T3.stat, tbin.pub, ord, tbin, pbdb.time, min.pub=10, plotit=makePlot))

##	corrected flux
pbdb.ord.flux <- apply(pbdb.ord.div,2,function(x) {
    raw.flux <- diff(c(0,x))
    return(raw.flux[raw.flux != 0])
})

##	sstat analysis
pbdb.sstat.ord.cor <- sstat.comp(pbdb.ord.flux,minN=10,plotit=makePlot)
