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

# make column for midpoint ma
pbdbNew$ma_mid <- (pbdbNew$max_ma + pbdbNew$min_ma) / 2

## use new instead of old
pbdb.dat <- pbdbNew

pbdbDat <- pbdbNew

## get rid of poor temporal resolution
pbdbDat <- pbdbDat[pbdbDat$tbin != '', ]

# get rid of bad taxonomy
pbdbDat <- pbdbDat[pbdbDat$family != '', ]
pbdbDat <- pbdbDat[pbdbDat$otu != '', ]

# get bin times
pbdbDat$mid_ma <- (pbdbDat$min_ma + pbdbDat$max_ma) / 2
pbdbTime <- sort(tapply(pbdbDat$mid_ma, pbdbDat$tbin, mean))
pbdbDat$tbin <- factor(pbdbDat$tbin, levels = names(pbdbTime))


# data.frame of publication, diversity and 3T stat
famTbinBias <- aggregate(list(div = pbdbDat$otu), list(fam = pbdbDat$family,
                                                       tbin = pbdbDat$tbin),
                         function(x) length(unique(x)))


## three timer stat

## matrix to determine three timers and part timers (sensu alroy 2008)
mt <- matrix(0, nrow = nlevels(pbdbDat$tbin), 
             ncol = nlevels(pbdbDat$tbin))
diag(mt) <- -10
mt[abs(row(mt) - col(mt)) == 1] <- 1

## loop through and compute three timers and part timers
timers <- lapply(split(pbdbDat$tbin, pbdbDat$otu), 
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

famTbinBias$T3Stat <- t3stat[match(famTbinBias$tbin, 
                                      levels(pbdbDat$tbin))]
famTbinBias$T3Div <- famTbinBias$div / famTbinBias$T3Stat

# record pubs per tbin
tbinPub <- tapply(pbdbDat$reference_no, pbdbDat$tbin,
                   function(x) length(unique(x)))
famTbinBias$tbinPub <- tbinPub[famTbinBias$tbin]


# calculate corrected diversity

pdf('ms/figSupp_divByPub.pdf', width = 4, height = 4)
pbdbFamDiv <- with(famTbinBias,
                   make3TPub(div, T3Stat, tbinPub, fam, tbin, pbdbTime, 
                             min.pub = 10, plotit = TRUE))
dev.off()

# corrected flux
pbdbFamFlux <- apply(pbdbFamDiv, 2, function(x) {
    flux <- diff(c(0, x)) 
    return(flux[flux != 0])
})

