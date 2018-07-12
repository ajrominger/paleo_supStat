#' @description function to produce a matrix of time by taxa with cells of corrected 
#' diversity
#' @param rawDiv
#' @param t3Stat
#' @param pub
#' @param taxaNames
#' @param tbin
#' @param tbinTime
#' @param minPub
#' @param plotit

make3tPub <- function(rawDiv, t3Stat, pub, taxaNames, tbin, tbinTime, minPub = 10, 
                      plotit = FALSE) {
	##	put data together so can be universally manipulated
	x <- data.frame(rawDiv = rawDiv, t3Stat = t3Stat, pub = pub, taxaNames = taxaNames, 
	                tbin = tbin)
	x$tbin <- as.character(x$tbin)
	x$taxaNames <- as.character(x$taxaNames)
	
	x <- x[!is.na(t3Stat) & pub >= minPub, ]
	
	tbinTime <- tbinTime[names(tbinTime) %in% x$tbin]
	
	##	3-timer correction
	t3Cor <- x$rawDiv/x$t3Stat
	
	##	publication correction
	logPub <- log(x$pub)
	pubMod <- lm(log(t3Cor) ~ logPub)
	pbdbPubMod <<- pubMod	# save regression to global env
	
	pub.resid <- exp(pubMod$residuals)
	
	##	plot so you can verify cuttoff etc.
	if(plotit) {
		plot(log(x$pub), log(t3Cor), 
			 xlab = 'log(Number of publications)',
			 ylab = 'log(3T-corrected number of genera)')
		
	    abline(pubMod, col = 'red')
	}
	
	tbinByTaxa <- socorro::tidy2mat(x$tbin, x$taxaNames, pub.resid)

	return(tbinByTaxa[names(sort(tbinTime, decreasing = TRUE)), ])
}
