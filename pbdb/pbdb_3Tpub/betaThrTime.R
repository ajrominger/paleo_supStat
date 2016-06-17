##  load sstat analysis
makePlot <- FALSE
source('~/Dropbox/Research/paleo_supStat/pbdb/pbdb_3Tpub/make_pbdb_3Tpub.R')

##  get first and last apperances
pbdb.ord.range <- do.call(rbind,with(ord.tbin.bias,tapply(tbin,ord,function(x) rev(range(pbdb.time[unique(x)])))))
pbdb.ord.range <- as.data.frame(pbdb.ord.range)
names(pbdb.ord.range) <- c("frst","last")

pbdb.ord.range$frst.bin <- names(pbdb.time[match(pbdb.ord.range[,1],pbdb.time)])
pbdb.ord.range$last.bin <- names(pbdb.time[match(pbdb.ord.range[,2],pbdb.time)])


##  combine time and beta for plotting ease
x <- data.frame(time=pbdb.ord.range$frst, beta=pbdb.sstat.ord.cor$beta[rownames(pbdb.ord.range)])
x <- x[!is.na(x$beta),]

source('~/R_functions/paleoAxis2.R')
source('~/R_functions/logAxis.R')

quartz(width=4,height=4)
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5,1,0))
plot(x,xlim=rev(range(x$time)),yaxt='n',log='y',
	xlab=expression('Time ('*10^6~'years)'),ylab=expression(beta))
abline(h=mean(pbdb.sstat.ord.cor$beta),lty=2,lwd=2)
logAxis(2)



