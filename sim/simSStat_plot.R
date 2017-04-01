library(socorro)
setwd('~/Dropbox/Research/paleo_supStat/sim')

sstatSim <- read.csv('simSStat.csv', as.is = TRUE)


plot(sstatSim$rho, sstatSim$var.shape, log = 'xy', axes = FALSE, type = 'n')
for(i in unique(sstatSim$S)) with(sstatSim[sstatSim$S == i & sstatSim$tmax == 100, ], 
                                  points(rho, var.shape, type = 'l', lty = 3))
for(i in unique(sstatSim$S)) with(sstatSim[sstatSim$S == i & sstatSim$tmax == 550, ], 
                                  points(rho, var.shape, type = 'l'))
logAxis(1)
logAxis(2)

plot(sstatSim$rho, sstatSim$var.rate, log = 'xy', axes = FALSE, type = 'n')
for(i in unique(sstatSim$S)) with(sstatSim[sstatSim$S == i & sstatSim$tmax == 100, ], 
                                  points(rho, var.rate, type = 'l', lty = 3))
for(i in unique(sstatSim$S)) with(sstatSim[sstatSim$S == i & sstatSim$tmax == 550, ], 
                                  points(rho, var.rate, type = 'l'))
logAxis(1)
logAxis(2)
