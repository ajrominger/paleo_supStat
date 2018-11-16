# **script to compare within clade fluctuation distributions at different taxonomic levels**

library(socorro)

# source needed functions
R.utils::sourceDirectory('R', modifiedOnly = FALSE)

load('data/pbdb_sstat_objects.RData')


# for each family, calculate aggregated eCDF and distribution of kurtosis values
famECDF <- lapply(sstatPBDBfam3TP$Px.sub, function(x) {
    simpECDF(scale(x)[, 1], complement = TRUE)
})
famECDF <- do.call(rbind, famECDF)
famKurt <- sapply(sstatPBDBfam3TP$Px.sub, kurt)

ordECDF <- lapply(sstatPBDBOrd$Px.sub, function(x) {
    simpECDF(scale(x)[, 1], complement = TRUE)
})
ordECDF <- do.call(rbind, ordECDF)
ordKurt <- sapply(sstatPBDBOrd$Px.sub, kurt)

clsECDF <- lapply(sstatPBDBCls$Px.sub, function(x) {
    simpECDF(scale(x)[, 1], complement = TRUE)
})
clsECDF <- do.call(rbind, clsECDF)
clsKurt <- sapply(sstatPBDBCls$Px.sub, kurt)

phyECDF <- lapply(sstatPBDBPhy$Px.sub, function(x) {
    simpECDF(scale(x)[, 1], complement = TRUE)
})
phyECDF <- do.call(rbind, phyECDF)
phyKurt <- sapply(sstatPBDBPhy$Px.sub, kurt)


#' @description function to plot theoretical and observed percentiles
#' @param x aggregated eCDF
#' @param ... additional plotting parameters

ppECDF <- function(x, ...) {
    alpha <- 0.75 / (1 + exp(0.0003 * (nrow(x) - 300))) # nicely scale transparency
    plot(pnorm(x[, 1], lower.tail = FALSE), x[, 2], pch = 16,
         col = gray(0, alpha = alpha), xlim = 0:1, ylim = 0:1, ...)
    
    abline(0, 1, col = 'red')
}


#' @description function to plot summary of kurtosis values distribution
#' @param x kurtosis values
#' @param ... additional plotting parameters

kurtInset <- function(x, ...) {
    allMean <- c(mean(famKurt), mean(ordKurt), mean(clsKurt), mean(phyKurt))
    allSD <- c(sd(famKurt), sd(ordKurt), sd(clsKurt), sd(phyKurt))
    
    plot(1, mean(x), pch = 16, 
         ylim = range(allMean - allSD, allMean + allSD) * c(1.5, 1.25), 
         xaxt = 'n', frame.plot = FALSE, yaxs = 'i', 
         ...)
    segments(x0 = 1, y0 = mean(x) - sd(x), y1 = mean(x) + sd(x))
}




# plot it
pdf('ms/figSupp_pkx_allTaxa.pdf', width = 9, height = 3)

split.screen(c(1, 4))

# first plots of the ECDF's
screen(1, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
ppECDF(famECDF)
mtext('Families', side = 3, line = 0.5)

screen(2, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
ppECDF(ordECDF, yaxt = 'n')
mtext('Orders', side = 3, line = 0.5)

screen(3, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
ppECDF(clsECDF, yaxt = 'n')
mtext('Classes', side = 3, line = 0.5)

screen(4, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
ppECDF(phyECDF, yaxt = 'n')
mtext('Phyla', side = 3, line = 0.5)


mtext('N(0, 1) percentiles', side = 1, outer = TRUE, line = 1.25)
mtext('Observed percentiles', side = 2, outer = TRUE, line = 1.25)

close.screen(all.screens = TRUE)


# now inset plots of kurtosis

start <- 1/4 + 0.01
swidth <- 1/32
increment <- 1/4 - 1/64
s <- split.screen(matrix(c(start + 0 * increment, start + 0 * increment + swidth, 0.25, 0.6, 
                           start + 1 * increment, start + 1 * increment + swidth, 0.25, 0.6, 
                           start + 2 * increment, start + 2 * increment + swidth, 0.25, 0.6, 
                           start + 3 * increment, start + 3 * increment + swidth, 0.25, 0.6),
                         ncol = 4, byrow = TRUE), erase = FALSE)

for(i in 1:4) {
    screen(s[i], new = FALSE)
    par(mar = rep(0, 4), mgp = c(1, 0.25, 0))
    kurtInset(switch(i, 
                     `1` = famKurt, 
                     `2` = ordKurt, 
                     `3` = clsKurt, 
                     `4` = phyKurt), 
              tcl = -0.1)
    mtext('Kurtosis', side = 2, line = 1.25)
}


close.screen(all.screens = TRUE)

dev.off()
