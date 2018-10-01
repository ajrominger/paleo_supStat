##	convenience function to produce a matrix of time by ord with cells of corrected diversity
make3TPub <- pbdb.3t.pub <- function(raw.div,t3.stat,pub,ord,tbin,tbin.time,min.pub=4,plotit=FALSE) {
	##	put data together so can be universally manipulated
	x <- data.frame(raw.div=raw.div,t3.stat=t3.stat,pub=pub,ord=ord,tbin=tbin)
	x$tbin <- as.character(x$tbin)
	x$ord <- as.character(x$ord)
	
	x <- x[!is.na(t3.stat) & pub >= min.pub,]
	
	tbin.time <- tbin.time[names(tbin.time) %in% x$tbin]
	
	##	3-timer correction
	t3.cor <- x$raw.div/x$t3.stat
	
	##	publication correction
	logPub <- log(x$pub)
	pub.lm <- lm(log(t3.cor)~logPub)
	pbdb.pub.lm <<- pub.lm	# save regression to global env
	
	pub.resid <- exp(pub.lm$residuals)
	
	##	plot so you can verify cuttoff etc.
	if(plotit) {
		plot(log(x$pub),log(t3.cor), 
			 xlab='log(Number of publications)',ylab='log(3T-corrected number of genera)')
		abline(pub.lm,col='red')
	}
	
	tbin.ord <- samp2site.spp(x$tbin,x$ord,pub.resid)

	return(tbin.ord[names(sort(tbin.time,decreasing=TRUE)),])
#	return(tbin.ord)
}
