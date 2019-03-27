library(socorro)

load('data/pbdb_sstat_objects.RData')
pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)

divFamRaw <- tapply(pbdbDat$otu, pbdbDat$family, function(x) length(unique(x)))

pdf('ms/figSupp_betaByRich.pdf', width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(divFamRaw[names(sstatPBDBfam3TP$beta)], sstatPBDBfam3TP$beta, log = 'xy', 
     xlab = 'Number of genera', ylab = expression(beta), 
     axes = FALSE, frame.plot = TRUE)
logAxis(1:2)
dev.off()


pdf('ms/figSupp_richGamma.pdf', width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(simpECDF(divFamRaw, complement = TRUE), log = 'x', xaxt = 'n', 
     xlab = 'Number of genera', ylab = 'Cumulative density')
logAxis(1)
param <- MASS::fitdistr(divFamRaw, 'gamma')
curve(pgamma(x, param$estimate[1], param$estimate[2], lower.tail = FALSE), col = 'red', add = TRUE)
legend('topright', legend = c('Observed', 'Best fit Gamma'), pch = c(1, NA), lty = c(NA, 1), 
       col = c('black', 'red'), bty = 'n')
dev.off()
