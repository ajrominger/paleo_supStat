########  methods for sstat objects  ########

# load library socorro if not already loaded
foo <- require('socorro')
if(!foo) devtools::install_github('ajrominger/socorro')

# likelihood function under sstat model
#' @param param the two parameters of the gamma distribution: shape and rate
#' @param dat the data with which to calculate the likelihood
sstatLLFun <- function(param, dat) {
	shape <- param[1]
	rate <- param[2]
	
	-sum(log(Px.gam(dat, shape, rate)))
}
	
# find the maximum likelihood estimate for sstat parameters
#' @param x an `sstat` object

sstatMLE <- function(x) {
	dat <- unlist(x$Px.sub)
	
	optim(c(0.55,0.17), sstatLLFun, method = 'BFGS', hessian = TRUE, dat = dat)
}


##	bootstrap likelihood
sstatBootMLE <- function(x, B = 1000) {
	dat <- x$Px.sub
	
	boots <- replicate(B, {
	    subDat <- sapply(dat, sample, size = 1)
	    
	    thisMLE <- try(optim(c(0.55, 0.17), sstatLLFun, method = 'BFGS', hessian = TRUE, dat = subDat), 
	                   silent = TRUE)
	    
	    if(class(this.mle) != 'try-error') {
	        if(thisMLE$convergence != 0) {
	            out <- rep(NA, 2)
	        } else {
	            out <- thisMLE$par
	        }
	    } else {
	        out <- rep(NA, 2)
	    }
	    
	    return(out)
	})
	
	# output
	o <- rbind(quantile(boots[1, ], c(0.025, 0.975), na.rm = TRUE),
	                  quantile(boots[2, ], c(0.025, 0.975), na.rm = TRUE))
	rownames(o) <- c('shape', 'rate')
	
	return(o)
}


# logLik for sstat object
#' @param x a `sstat` object
#' @param fitted logical, whether the super-stats prediction was arrived at by 
#' fitting with ML, or by integrating over empirical p_k(x|beta) and f(beta) 
#' distributions

logLik.sstat <- function(x, fitted = FALSE) {
	dat <- unlist(x$Px.sub)
	
	lik <- sum(log(x$Px(dat)))
	
	if(fitted) {
		attr(lik, 'df') <- 2
	} else {
		attr(lik, 'df') <- 0
	}
	
	class(lik) <- 'logLik'
	
	return(lik)
}


# plot for sstat object
#' @param x an `sstat` object
#' @param sstatCol the plotting color for super-stats prediction
#' @param normCol the plotting color for the Gaussian prediction
#' @param showNorm logical, whether or not to show the Gaussian prediciton
#' @param addLegend logical, whether or not to add a legend
#' @param ... further arguments passed to `plot`

plot.sstat <- function(x, sstatCol = 'red', normCol = 'blue', showNorm = TRUE, 
                       addLegend = TRUE, ...) {
	dots <- list(...)
    thisCDF <- socorro::simpECDF(abs(unlist(x$Px.sub)), complement = TRUE)
	
	plot(thisCDF, log='xy', #col='gray',pch=16,
		 ...)
	
	if(dots$xaxt != 'n' | !('xaxt' %in% names(dots))) socorro::logAxis(1)
	if(dots$yaxt != 'n' | !('yaxt' %in% names(dots))) socorro::logAxis(1)
	
	PPx <- x$PPx
	curve(PPx(x, comp = TRUE), col = sstatCol, lwd = 2, add = TRUE)
	
	if(showNorm) {
		thisSD <- sd(unlist(x$Px.sub))
		curve(2 * pnorm(x, 0, thisSD, lower.tail = FALSE), col = normCol, lwd = 2, add = TRUE)
	}
	
	if(addLegend) {
		leg <- c('Data', 'Superstatistics')
		col <- c(par('fg'), sstatCol)
		pch <- c(1, NA)
		pt.cex <- c(1, NA)
		lty <- c(NA, 1)
		lwd <- c(NA, 2)
		
		if('panel.first' %in% names(list(...))) {
			leg <- c(leg, 'Superstatistics likelihood CI')
			col <- c(col, socorro::colAlpha(sstatCol, 0.5))
			pt.cex <- c(pt.cex, 2)
			lwd <- c(lwd, NA)
			lty <- c(lty, NA)
			pch <- c(pch, 15)
		}
		
		if(showNorm) {
			leg <- c(leg, 'ML normal')
			col <- c(col, normCol)
			pt.cex <- c(pt.cex, NA)
			lwd <- c(lwd, 2)
			lty <- c(lty, 1)
			pch <- c(pch, NA)
			
			# if('panel.first' %in% names(list(...))) {
			# 	leg <- c(leg,'Normal likelihood CI')
			# 	col <- c(col,hsv(0.6, s=0.5))
			# 	pt.lwd <- c(pt.lwd,1)
			# 	pt.cex <- c(pt.cex,2)
			# 	lwd <- c(lwd,NA)
			# 	pch <- c(pch,15)
			# }
		}
		
		extracex <- ifelse('cex' %in% names(list(...)), list(...)$cex, 1)
		legend('bottomleft', legend = leg, col = col, pch = pch, 
		       pt.cex = pt.cex * extracex, lwd = lwd, lty = lty, bty='n')
	}
}


# function to add confidence interval polygon from ML analysis
#' @param ci a matrix of ....
#' @param fun ....
#' @param ... additional parameters passed to `polygon`

mlePoly <- function(ci, fun, ...) {
	n <- 100
	
    x <- 10^c(seq(par('usr')[1], par('usr')[2], length.out = n), 
	          seq(par('usr')[2], par('usr')[1], length.out = n))
	
    y <- c(fun(x[1:n], ci[1, 1], ci[2, 2]), fun(x[(1:n) + n], ci[1, 2], ci[2, 1]))
	y[y < 10^par('usr')[3]] <- 10^par('usr')[3]
	
	polygon(x = x, y = y, ...)
}

