library(socorro)
setwd('~/Dropbox/Research/paleo_supStat/sim')

## helper function for plotting *s*uper *s*at *s*imulation
sssPlot <- function(x, y, log, tt, ttlty) {
    spp <- unique(sstatSim$S)
    nspp <- length(spp)
    sppCol <- hsv(seq(0.45, 0.7, length.out = nspp), 
                  seq(0.5, 1, length.out = nspp), 
                  seq(1, 0.75, length.out = nspp))
    
    plot(sstatSim[, x], sstatSim[, y], log = log, type = 'n', 
         xaxt = ifelse(grepl('x', log), 'n', 's'), 
         yaxt = ifelse(grepl('y', log), 'n', 's'), 
         ...)
    
    for(j in 1:length(tt)) {
        for(i in unique(sstatSim$S)) 
            with(sstatSim[sstatSim$S == i & sstatSim$tmax == tt[j], ], 
                 points(get(x), get(y), type = 'l', lty = ttlty[j], lwd = 2,
                        col = quantCol(i, sppCol, xlim = range(spp))))
    }
    
    if(grepl('x', log)) logAxis(1)
    if(grepl('y', log)) logAxis(2)
}

sstatSim <- read.csv('simSStat.csv', as.is = TRUE)


sssPlot('rho', 'var.shape', log = 'xy', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'var.rate', log = 'xy', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'mean.mean', log = '', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'mean.sd', log = '', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'minS', log = 'xy', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'maxS', log = 'xy', tt = c(100, 550), ttlty = c(3, 1))

sssPlot('rho', 'endT', log = 'y', tt = 550, ttlty = 1)
