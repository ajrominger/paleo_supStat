# **script to compare within clade volatility distributions at different taxonomic levels**

library(socorro)

# source needed functions
R.utils::sourceDirectory('R', modifiedOnly = FALSE)

load('data/pbdb_sstat_objects.RData')


#' @description function to plot f(beta) distribution
#' @param obj the sstat object
#' @param thrCol color for plotting of theoretical curve

fbetaPlot <- function(obj, thrCol = 'red', ...) {
    betaCDF <- simpECDF(obj$beta, complement = TRUE)
    plot(betaCDF, ylim = c(0, 1.025),
         log = 'x', xaxt = 'n', yaxt = 'n',
         xlab = expression(beta), ...)
    
    logAxis(1, expLab = TRUE)
    
    curve(pgamma(x, obj$gam.par[1], obj$gam.par[2], 
                 lower.tail = FALSE), 
          col = thrCol, lwd = 2, add = TRUE)
    
}


# the plotting
pdf('ms/figSupp_fbeta_allTaxa.pdf', width = 9, height = 3)

split.screen(c(1, 4))

screen(1, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
fbetaPlot(sstatPBDBfam3TP)
axis(2)
mtext('Families', side = 3, line = 0.5)

screen(2, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
fbetaPlot(sstatPBDBOrd)
mtext('Oders', side = 3, line = 0.5)

screen(3, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
fbetaPlot(sstatPBDBCls)
mtext('Classes', side = 3, line = 0.5)

screen(4, new = FALSE)
par(mar = c(0.3, 0.3, 1.5, 0.3), oma = c(2.5, 2.5, 0, 0), 
    mgp = c(2, 0.5, 0))
fbetaPlot(sstatPBDBPhy)
mtext('Phyla', side = 3, line = 0.5)

mtext(expression(beta), side = 1, outer = TRUE, line = 1.25)
mtext('Cumulative density', side = 2, outer = TRUE, line = 1.25)

close.screen(all.screens = TRUE)
dev.off()
