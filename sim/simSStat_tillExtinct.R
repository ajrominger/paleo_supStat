library(parallel)
library(MASS)
library(socorro)
setwd('~/Dropbox/Research/paleo_supStat/sim')

## function to time average origination and extinction
tAvg <- function(tt, oe, tbin, tmax) {
    oe <- oe[tt < tmax]
    tt <- tt[tt < tmax]
    tbins <- seq(0, tmax, by = tbin)
    
    ## remove empty tbins
    tbins <- tbins[1:(max(which(tbins <= max(tt))) + 1)]
    
    ## divide the sequence by tbins and average
    tcut <- cut(tt, tbins)
    
    avg <- sapply(split(oe, tcut), function(x) c(sum(x == 1), sum(x ==-1)))
    
    return(avg)
}

## control parameters on simulation
nrate <- 10
tbin <- 10

allPar <- expand.grid(exp(seq(log(0.001), log(3), length.out = nrate)),
                      seq(10, 100, length.out = ceiling(nrate/3)))
allRho <- allPar[, 1]
allS <- allPar[, 2]
tmax <- 1e+03

sstatSim <- mclapply(1:length(allRho), 
                     mc.cores = 6,
                     FUN = function(i) {
    print(i)
    
    ## set simulation parameters
    la <- mu <- allRho[i]
    S <- allS[i]
    
    ## calculate an over-estimate of number of timesteps till hitting tmax
    nt <- 10 * tmax / (1/(S*(la+mu)))
    
    ## simulate multiple realizations
    out <- replicate(500, {
        oe <- c(0, sample(c(-1, 1), nt-1, replace = TRUE))
        St <- S + cumsum(oe)
        
        ## truncate the realization to extinction event
        if(any(St <= 0)) {
            oe <- oe[1:(min(which(St <= 0)))]
            St <- St[1:(min(which(St <= 0)))]
        }
        
        ## simulate time, always starting at time 0
        tt <- cumsum(c(0, rexp(length(St)-1, St[-length(St)] * (la + mu))))
        
        ## if 2 or less time bins, return NAs
        if(ceiling(max(tt) / tbin) < 3) {
            out <- c(mean = NA, var = NA)
        } else {
            ## calculate the time average of orig and extinction
            oeAvg <- tAvg(tt, oe, tbin, tmax)
            
            ## return mean of delta St and its var
            out <- c(mean = mean(oeAvg[1, -1] - oeAvg[2, -ncol(oeAvg)]), 
                     var = var(oeAvg[1, ]) + var(oeAvg[2, ]))
        }
        
        ## add to output the range of S, whether the process went extinct (`ext = 1` 
        ## coresponds to extinct) and what time the process ends at
        out <- c(out, Smin = min(St[St > 0]), Smax = max(St), 
                 ext = as.numeric(St[length(St)] == 0), 
                 endT = max(tt[tt <= tmax]))
        return(out)
    })
    
    gpar <- try(fitdistr(1/out[2, out[2, ] > 0 & !is.na(out[2, ]) & out[5, ] == 1], 'gamma'), 
                silent = TRUE)
    if(class(gpar) == 'try-error') gpar <- list(estimate = c(shape = NA, rate = NA))
    
    mpar <- c(mean = mean(out[1, !is.na(out[1, ])]), sd = sd(out[1, !is.na(out[1, ])]))
    ext <- mean(out[5, ] == 1)
    endT <- mean(out[6, ])
    minS <- mean(out[3, ])
    maxS <- mean(out[4, ])
    
    return(c(mean = mpar, var = gpar$estimate, minS = minS, maxS = maxS, 
             ext = ext, endT = endT))
    
})

sstatSim <- do.call(rbind, sstatSim)

## combine with sim params and write out
sstatSim <- cbind(rho = allRho, S = allS, tmax = tmax, sstatSim)
write.csv(sstatSim, file = 'simSStat_tillExtinct.csv', row.names = FALSE)
