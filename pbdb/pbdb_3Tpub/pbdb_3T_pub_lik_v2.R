oldwd <- setwd('~/Dropbox/Research/paleo_supStat')

##	needed functions and data
source("code/sstat_methods.R")

if(!exists("pbdb.sstat.ord.cor")) {
	cat('making sstat objects','\n')
	makePlot <- FALSE
	source("pbdb/pbdb_3Tpub/make_pbdb_3Tpub.R")
}

##	ml confidence bootstrap interval
if(!exists("pbdb.ord.mleCI")) {
	cat('constructing ML CI','\n')
	pbdb.ord.mleCI <- boot.mle.sstat(pbdb.sstat.ord.cor,B=1000,use.all=FALSE)
}

##	plot it
plot(pbdb.sstat.ord.cor,xlab="|Fluctuations|",ylab="Cumulative density", cex=2, 
	 panel.first={
	 	mle.poly(pbdb.ord.mleCI$norm,function(x,mean,sd) {2*pnorm(x,mean,sd,lower.tail=FALSE)},
	 			 col=hsv(0.6, 0.5),border=NA)
	 	mle.poly(pbdb.ord.mleCI$sstat,PPx.gam,
	 			 col=hsv(s=0.5),border=NA)
})

##  return to original wd
setwd(oldwd)