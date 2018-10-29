#' @description function to produce a matrix of time by taxa with cells of corrected diversity
#' @param rawDiv the raw diversity of each taxon in each time interval
#' @param t3stat the 3 timer stat for each diversity record
#' @param pub the number of publications associated with each diversity record
#' @param taxa the taxon names for each diversity record
#' @param tbin the time interval of each diversity record
#' @param tbinTime times associated with each `tbin`
#' @param minPub minimum number of publications for inclusion in regression analysis
#' @param plotit logical, should plot of taxon richness versus number of publications be made
#' @return a matrix with rows corresponding to time intervals and columns to the given taxa
#' each cell in the matrix represents corrected taxon richness


make3TPub <- function(rawDiv,  t3stat,  pub,  taxa,  tbin,  tbinTime,  
                      minPub = 10,  plotit = FALSE) {
	# put data together so can be universally manipulated
	x <- data.frame(rawDiv = rawDiv, t3stat = t3stat, pub = pub, taxa = taxa, tbin = tbin)
	x$tbin <- as.character(x$tbin)
	x$taxa <- as.character(x$taxa)
	
	x <- x[!is.na(t3stat) & pub >=  minPub, ]
	
	tbinTime <- tbinTime[names(tbinTime) %in% x$tbin]
	
	# 3-timer correction
	t3cor <- x$rawDiv/x$t3stat
	
	# publication correction
	logPub <- log(x$pub)
	pubLM <- lm(log(t3cor)~logPub)
	pbdbPubLM <<- pubLM	# save regression to global env
	
	pubResid <- exp(pubLM$residuals)
	
	# plot so you can verify cuttoff etc.
	if(plotit) {
		plot(log(x$pub), log(t3cor),  
			 xlab = 'log(Number of publications)', 
			 ylab = 'log(3T-corrected number of genera)')
		abline(pubLM, col = 'red')
	}
	
	tbinTaxa <- socorro::tidy2mat(x$tbin, x$taxa, pubResid)

	return(tbinTaxa[names(sort(tbinTime, decreasing = TRUE)), ])
}
