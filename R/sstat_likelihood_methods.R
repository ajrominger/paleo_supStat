#' @description gives the log likelihood function under sstat model
sstatLL <- function(par, dat) {
    -sum(log(Px.gam(dat,par[1], par[2])))
}

##	find the maximum likelihood estimate for sstat parameters
MLE.sstat <- function(x, useAll = FALSE) {
    if(useAll) {
        theseDat <- unlist(x$Px.raw)
    } else {
        theseDat <- unlist(x$Px.sub)
    }
    
    optim(c(0.55, 0.17), sstatLL, method = 'BFGS', hessian = TRUE, dat = theseDat)
}

#' @description bootstrap likelihood for super stats model
#' @param x the `sstat` object
#' @param B the number of boostrap replicates
#' @param useAll logical, whether all orders, or only those with the minimum number of occurences as specified
#' in `make3TPub` argument `minPub` should be used

bootMLE.sstat <- boot.mle.sstat <- function(x, B = 1000, useAll = FALSE) {
    if(useAll) {
        theseDat <- x$Px.raw
    } else {
        theseDat <- x$Px.sub
    }
    
    boots <- replicate(B, {
        # browser()
        subDat <- sapply(theseDat, sample, size = 1)
        thisMLE <- try(MLE.sstat(x, useAll), silent = TRUE)
        # thisMLE <- try(optim(c(0.55, 0.17), sstat.lik, method = 'BFGS', hessian = TRUE, dat = subDat), 
        #                silent = TRUE)
        
        if(class(thisMLE) !=  'try-error') {
            if(thisMLE$convergence !=  0) {
                out <- rep(NA, 2)
            } else {
                out <- thisMLE$par
            }
        } else {
            out <- rep(NA, 2)
        }
        out
    })
    
    sstatOut <- rbind(quantile(boots[1, ], c(0.025, 0.975), na.rm = TRUE), 
                      quantile(boots[2, ], c(0.025, 0.975), na.rm = TRUE))
    rownames(sstatOut) <- c('shape', 'rate')
    
    return(list(sstat = sstatOut))
}
