source('~/R_functions/my_ecdf.R')
source('~/R_functions/normLS.R')
source('~/R_functions/gammaLS.R')
source('~/R_functions/logPlot.R')
source('~/Dropbox/Research/paleo_supStat/code/Px_gam.R')

# function to make `sstat` object
#' @param dat the data to be used, organized as a list of taxa, each component of the list
#' being a vector of that taxon's fluctuations through time
#' @param minN the minimum number of observations to use, taxa with fewer observations are ignored
#' @param plotit logical, should a plot be made
#' @param ... arguments passed to plotting function if plot is to be made

sstatMake <- function(dat, minN = 15, plotit = TRUE, ...) {
	
    these2use <- sapply(dat, length) >= minN
	p2use <- dat[these2use]
	
	cat('computing Gaussian fit for p_k(x|beta) \n')
	pkPar <- sapply(p2use, function(x) unlist(normLS(x)[c('par', 'value')]))
	pkPar <- t(pkPar)
	colnames(pkPar) <- c('mu', 'sig', 'ss')
	
	cat('re-centering \n')
	for(i in 1:length(p2use)) {
		p2use[[i]] <- p2use[[i]] - pkPar[i, 'mu']
	}
	
	cat('computing f(beta) \n')
	fPar <- gammaLS(1 / (pkPar[, 'sig'])^2)$par
	fuentPar <- c(n = 2 * fPar[1], b0 = fPar[1] * fPar[2])
	
	cat('computing P(x) \n')
	thisPx <- function(x) Px.gam(x, fPar[1], fPar[2])
	thisPPx <- function(x, comp = TRUE) PPx.gam(x, fPar[1], fPar[2], comp)
	
	
	# object to return
	out <- list(gamPar = fPar, sspar = fuentPar, beta = 1 / (pkPar[, 'sig'])^2, 
	            sumSq = pkPar[, 'ss'], minN = minN, Pxsub = p2use, Px = thisPx, 
	            PPx = thisPPx)
	
	class(out) <- 'sstat'
	
	if(plotit) {
		plot.sstat(out, ...)
	}
	
	return(out)
}

