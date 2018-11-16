# **a script to compare our 3TPub curve to other estimates of richness through the Phanerozoic**

# package with diversity dynamics subsampling functions
library(divDyn)

# package for plotting
library(socorro)

# load and prep data
pbdbFamDiv <- read.csv('data/pbdb_3TPub-corrected.csv', row.names = 1)
pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)
tbin <- read.csv('data/tbinsMid.csv', as.is = TRUE)
tbin$tbin <- factor(tbin$tbin, levels = tbin$tbin)
pbdbDat$tbin <- factor(pbdbDat$tbin, levels = levels(tbin$tbin))
pbdbDat$tbinNum <- as.integer(pbdbDat$tbin)

pbdbDatUnique <- pbdbDat[!duplicated(paste0(pbdbDat$collection_no, pbdbDat$otu)), ]

# subsampled richness
pbdbCR <- subsample(pbdbDatUnique, bin = 'tbinNum', tax = 'otu', iter = 10, q = 55, 
                    type = 'cr', unit = 'reference_no')
pbdbSQS <- subsample(pbdbDatUnique, bin = 'tbinNum', tax = 'otu', iter = 50, q = 0.75, 
                     ref = 'reference_no', type = 'sqs', singleton = 'ref')

# our new richness estimate
pbdbT3Pub <- rowSums(pbdbFamDiv)


# plot fluctuations to see that they're comprable
pdf('ms/figSupp_divEstComp.pdf', width = 7.5, height = 4)
layout(matrix(1:2, nrow = 1))

par(mar = c(4.5, 2.5, 0, 0.5) + 0.5, mgp = c(2, 0.75, 0))
plot(1, xlim = c(540, 0), ylim = c(-500, 500), type = 'n', xaxt = 'n', 
     xlab = '', ylab = 'Richness', xaxs = 'i')
paleoAxis(1)
mtext('Millions of years ago', side = 1, line = 3.5)

lines(tbin$ma_mid[-1], diff(pbdbCR$divCSIB), col = 'black', lwd = 2)
lines(tbin$ma_mid[-1], diff(pbdbSQS$divCSIB), col = 'blue', lwd = 2)
lines(tbin$ma_mid[-c(1:2, nrow(tbin))], diff(pbdbT3Pub), col = 'red', lwd = 2)


par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(simpECDF(abs(diff(pbdbT3Pub)), complement = TRUE), col = 'red', log = 'xy', 
     type = 'l', lwd = 2, xlim = c(1, 500),
     panel.first = {
         lines(simpECDF(abs(diff(pbdbCR$divCSIB)), complement = TRUE), 
               col = 'black', lwd = 2)
         lines(simpECDF(abs(diff(pbdbSQS$divCSIB)), complement = TRUE), 
               col = 'blue', lwd = 2)
     }, 
     axes = FALSE, frame.plot = TRUE, 
     xlab = '|Fluctuations|', ylab = 'Cumulative density')
logAxis(1:2)

dev.off()
