pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)
famDiv <- read.csv('data/pbdb_3TPub-corrected.csv', row.names = 1)



divFamRaw <- tapply(pbdbDat$otu, pbdbDat$family, function(x) length(unique(x)))
divFamCor <- colSums(famDiv)
foo <- as.matrix(famDiv)
foo[foo == 0] <- NA

divFamMean <- colMeans(foo, na.rm = TRUE)

plot(1/sstatPBDBfam3TP$beta, divFamCor[names(sstatPBDBfam3TP$beta)], log = 'xy')
plot(1/sstatPBDBfam3TP$beta, divFamMean[names(sstatPBDBfam3TP$beta)], log = 'xy')
