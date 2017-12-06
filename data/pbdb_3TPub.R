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
occsRaw <- txg

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

## remove taxa that have occ only in first and last bins
colnames(occsFill) <- colnames(occsRaw)
occsFill <- occsFill[, colnames(occsFill) %in% unique(pbdb$otu)]
occsRaw <- occsRaw[, colnames(occsRaw) %in% unique(pbdb$otu)]

## correct raw occs by 3 timer
occs3T <- occsRaw / presProb
dim(occsRaw)
length(presProb)

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

## plot of 3timer corrected genus richness through time
par(xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, rowSums(occsRaw) / presProb, type = 'l', xlim = c(540, 0),
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)

par(xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, rowSums(occs3T), type = 'l', xlim = c(540, 0),
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)


## correct for number of publications
logD <- log(rowSums(occs3T))
logP <- log(tapply(pbdb$reference_no, pbdb$tbin, function(x) length(unique(x))))

pubMod <- lm(logD ~ logP)
pubCorFact <- as.numeric(pubMod$coefficients[1] * exp(logP) ^ pubMod$coefficients[2])
occs3TPub <- occs3T / pubCorFact

par(mar = c(2.5, 2.5, 0, 0) + 0.5, xpd = FALSE, mgp = c(1.75, 0.5, 0))
plot(logP, logD)
abline(pubMod, col = 'red')

par(mar = c(4, 2.5, 0, 0) + 0.5, xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, logD, type = 'l', xlim = c(540, 0), ylim = range(logP, logD),
     xaxt = 'n', xaxs = 'i')
lines(tbins$ma_mid, logP, col = 'red')
paleoAxis(1)

par(mar = c(4, 2.5, 0, 0) + 0.5, xpd = NA, mgp = c(1.5, 0.5, 0))
plot(tbins$ma_mid, rowSums(occs3TPub), type = 'l', xlim = c(540, 0),
     xaxt = 'n', xaxs = 'i')
paleoAxis(1)

## corect for number of pubs, but at order level
occsNames2Ord <- pbdb$order[match(colnames(occs3T), pbdb$otu)]
occs3T <- occs3T[, occsNames2Ord != '']
occsNames2Ord <- occsNames2Ord[occsNames2Ord != '']

occs3TbyOrd <- split(as.data.frame(t(occs3T)), occsNames2Ord)
ordDiv <- lapply(occs3TbyOrd, colSums)
logPOrd <- rep(logP, length(ordDiv))
logOrdDiv <- log(unlist(ordDiv))

plot(logPOrd, logOrdDiv)
