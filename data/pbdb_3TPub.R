library(socorro)

setwd('~/Dropbox/Research/paleo_supStat/data')

## read in pbdb data and time bins
tbins <- read.csv('timebins.csv', as.is = TRUE)
tbins$ma_mid <- (tbins$ma_min + tbins$ma_max) / 2
pbdb <- read.csv('pbdb_data.csv', as.is = TRUE)
pbdb$tbin <- factor(pbdb$tbin, levels = tbins$tbin)

## matrix of time bin by genus
txg <- tidy2mat(pbdb$tbin, pbdb$otu, rep(1, nrow(pbdb)))

## range throughs and gaps
gapRangeThrough <- sapply(1:ncol(txg), function(i) {
    x <- txg[, i]
    occs <- which(x > 0)
    nots <- which(x == 0)
    holes <- nots[nots > min(occs) & nots < max(occs)]
    
    rangeThroughs <- gaps <- occsFill <- rep(0, nrow(tbins))
    
    ## only record range throughs and gaps if there could have been a gap
    if(length(occs) > 2 | length(holes) > 0) {
        rangeThroughs[min(occs):max(occs)] <- 1
        gaps[holes] <- 1
    }
    
    ## fill in gaps
    occsFill[min(occs):max(occs)] <- 1
    
    return(c(gaps, rangeThroughs, occsFill))
})

## tbin by genus matrix will filled-in durations
occsFill <- gapRangeThrough[(1:nrow(tbins)) + 2 * nrow(tbins), ]
occsRaw <- gapRangeThrough[(1:nrow(tbins)) + nrow(tbins), ]

## probability of preservation
gapRangeThrough <- rowSums(gapRangeThrough)
presProb <- gapRangeThrough[1:nrow(tbins)] / gapRangeThrough[(1:nrow(tbins)) + nrow(tbins)]

## remove first and last time bins
goodTimes <- 2:(nrow(tbins) - 1)
occsFill <- occsFill[goodTimes, ]
occsRaw <- occsRaw[goodTimes, ]
presProb <- presProb[goodTimes]
tbins <- tbins[goodTimes, ]
pbdb <- pbdb[pbdb$tbin %in% tbins$tbin, ]
pbdb$tbin <- factor(as.character(pbdb$tbin), levels = tbins$tbin)

## plot of preservation probability through time
par(xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, presProb, type = 'l', xlim = c(540, 0), ylim = 0:1,
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)

## plot of raw genus richness through time
par(xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, rowSums(occsFill), type = 'l', xlim = c(540, 0),
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)

## plot of 3timmer corrected genus richness through time
par(xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, rowSums(occsRaw) / presProb, type = 'l', xlim = c(540, 0),
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)

'Jurassic 6' %in% pbdb$tbin

tapply(pbdb$reference_no, pbdb$tbin, function(x) length(unique(x)))
