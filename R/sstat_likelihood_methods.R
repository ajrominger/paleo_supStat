#' @description gives the log likelihood function under sstat model
#' @param par the parameter values
#' @param dat the data

sstatLL <- function(par, dat) {
    -sum(log(Px.gam(dat,par[1], par[2])))
}


#' @description finds the maximum likelihood estimate of the superstats model
#' @param dat the data to fit

sstatMLE <- function(dat) {
    optim(c(0.55, 0.17), sstatLL, method = 'BFGS', hessian = TRUE, 
          dat = dat)
}


#' @description bootstrap likelihood for super stats model
#' @param x the `sstat` object
#' @param B the number of boostrap replicates
#' @param useAll logical, whether all orders, or only those with the minimum number of 
#' occurences as specified
#' in `make3TPub` argument `minPub` should be used

bootMLE.sstat <- function(x, B = 1000, useAll = FALSE) {
    if(useAll) {
        theseDat <- x$Px.raw
    } else {
        theseDat <- x$Px.sub
    }
    
    boots <- replicate(B, {
        subDat <- sapply(theseDat, sample, size = 1)
        
        thisMLE <- try(sstatMLE(subDat), silent = TRUE)
        
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


#' @description logLik for sstat class
#' @param x the `sstat` object
#' @param fitted logical, was the model fitted by max likelihood or computed from first 
#' principles
#' @param useAll logical, should all data be used, or only those taxa that have greater 
#' than `minN` occurrences
#' as specified in `sstatComp`

logLik.sstat <- function(x, fitted = TRUE, useAll = FALSE) {
    if(useAll) {
        theseDat <- unlist(x$Px.raw)
    } else {
        theseDat <- unlist(x$Px.sub)
    }
    
    lik <- sum(log(x$Px(theseDat)))
    
    if(fitted) {
        attr(lik, 'df') <- 2
    } else {
        attr(lik, 'df') <- 0
    }
    
    class(lik) <- 'logLik'
    
    return(lik)
}
