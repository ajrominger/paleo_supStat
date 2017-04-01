library(socorro)
library(parallel)
library(MASS)

## function to time average origination and extinction
tAvg <- function(tt, oe, tbin, tmax) {
    oe <- oe[tt < tmax]
    tt <- tt[tt < tmax]
    tbins <- seq(0, tmax, by = tbin)
    
    ## remove empty tbins
    tbins <- tbins[1:(max(which(tbins <= max(tt))) + 1)]
    
    ## divide the sequence by tbins and average
    tcut <- try(cut(tt, tbins))
    if(class(tcut) == 'try-error') browser()
    
    avg <- sapply(split(oe, tcut), function(x) c(sum(x == 1), sum(x ==-1)))
    
    return(avg)
}

## control parameters on simulation
ntaxa <- 10
tmax <- 550
tbin <- 10

allRho <- exp(seq(log(0.001), log(2), length.out = ntaxa))
allS <- rep(100, ntaxa)
allTmax <- rep(550, ntaxa)

gammaPar <- lapply(1:ntaxa, 
                     # mc.cores = 6,
                     FUN = function(i) {
    ## set simulation parameters
    la <- mu<- allRho[i]
    S <- allS[i]
    
    ## calculate an over-estimate of number of timesteps till hitting tmax
    nt <- 10 * tmax / (1/(S*(la+mu)))
    
    # browser()
    
    ## simulate multiple realizations
    out <- replicate(50, {
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
        if(ceiling(tt / tbin) < 3) {
            out <- c(mean = NA, var = NA)
        } else {
            ## calculate the time average of orig and extinction
            oeAvg <- tAvg(tt, oe, tbin, tmax)
            
            ## return mean of delta St and its var
            out <- c(mean = mean(oeAvg[1, -1] - oeAvg[2, -ncol(oeAvg)]), 
                     var = var(oeAvg[1, ]) + var(oeAvg[2, ]))
        }
        
        ## add to output whether the process went extinct (`ext = 1` coresponds to extinct)
        ## and what time the process ends at
        out <- c(out, ext = as.numeric(St[length(St)] == 0), endT = max(tt))
        return(out)
    })
    
    gpar <- fitdistr(1/out[2, ], 'gamma')
    mpar <- c(mean = mean(out[1, ]), sd = sd(out[1, ]))
    ext <- mean(out[3, ] == 1)
    
    return(c(mean = mpar, var = gpar$estimate, ext = ext))
})

gammaPar <- do.call(rbind, gammaPar)

gammaPar <- cbind(rho = allRho, S = allS, tmax = allTmax, gammaPar)



plot(allRho, gammaPar[, 1], log = 'xy')
plot(allRho, gammaPar[, 2], log = 'xy')
