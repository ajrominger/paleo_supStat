##	need to figure out what this depends on...

library(fields)
library(classInt)


pbdb.ord.range <- do.call(rbind,with(ord.tbin.bias,tapply(tbin,ord,function(x) rev(range(pbdb.time[unique(x)])))))
pbdb.ord.range <- as.data.frame(pbdb.ord.range)
names(pbdb.ord.range) <- c("frst","last")

pbdb.ord.range$frst.bin <- names(pbdb.time[match(pbdb.ord.range[,1],pbdb.time)])
pbdb.ord.range$last.bin <- names(pbdb.time[match(pbdb.ord.range[,2],pbdb.time)])

pbdb.time.col <- tim.colors(length(pbdb.time))
names(pbdb.time.col) <- names(sort(pbdb.time))
pbdb.time.col <- pbdb.time.col[names(pbdb.time)]

pbdb.sstat.ord.cor$beta

paleoPlot(pbdb.ord.range$frst,log(pbdb.sstat.ord.cor$beta[rownames(pbdb.ord.range)]))

##	function to plot beta through time
source("~/R_functions/paleoPlot.R")
source("~/R_functions/logAxis.R")
plot.betTime <- function(b,tm,col,frst=TRUE,...) {
	these.time <- tm[names(b),]
	
	if(missing(col)) col <- rep(ifelse(frst,"blue","red"),length(b))
	
	if(!frst) these.col[these.time[,2] == min(these.time[,2])] <- NA
	
	plot(these.time[,ifelse(frst,1,2)],b,
			  # log="y",
			  # xlim=rev(range(these.time[,1:2])),
			  xlim=c(550,0),
			  col=col,
			  xlab=ifelse(frst,"Origination time","Extinction time"),...)
	# text(these.time[,ifelse(frst,1,2)],b,labels=names(b))
	# abline(h=mean(b))
}

palette(tim.colors(length(unique(these.cls))))
plot.betTime(pbdb.sstat.ord.cor$beta,pbdb.ord.range,ylab=expression(beta),
			 col=as.factor(these.cls),
			 log="y",yaxt="n")
logAxis(2)
text(c(240,260),c(0.004,0.02),labels=c("Ammonoidae","Scleractinia"),pos=4)

par(oma=rep(0.8,4),mar=rep(0,4),mfrow=c(3,5))

for(i in 1:length(which(table(these.cls)>1))) {
	incld.ths <- which(these.cls == names(which(table(these.cls)>1))[i] & !(names(pbdb.sstat.ord.cor$beta) %in% c("Scleractinia","Ammonoidea")))
	plot.betTime(pbdb.sstat.ord.cor$beta[incld.ths],pbdb.ord.range,
				 ylab=expression(beta),
				 log="y",yaxt="n",xaxt="n",ylim=c(0.004,40))
	abline(h=mean(pbdb.sstat.ord.cor$beta),lty=2)
	legend("bottomright",legend=names(which(table(these.cls)>1))[i],bty="n")
}



##	which classes have gone extinct
tapply(pbdb.ord.range[names(pbdb.sstat.ord.cor$beta),"last"],these.cls,min)



plot.betTime(pbdb.sstat.ord.cor$beta,pbdb.ord.range,ylab=expression(beta),frst=FALSE,yaxt="n")
logAxis(2)
text(70,0.004,labels="Ammonoidae",pos=2)





















plot.betRange <- function(b,tm,...) {
	these.time <- tm[names(b),]
	
	these.col <- pbdb.time.col[as.character(these.time[,3])]
	
	these.col[these.time[,2] == min(these.time[,2])] <- "transparent"
	
	plot(these.time[,1]-these.time[,2],b,
			  #col=these.col,
			  ...)
	# abline(h=mean(b))
	# legend("bottomright",legend=min(these.time[,1]),pch=1,col=hsv(min(these.time[,1]/max(these.time[,1]))))
}

with(pbdb.sstat.ord.cor,
	 plot.betRange(beta,pbdb.ord.range,
				   xlab="Order lifetime",ylab=expression(beta),
				   yaxt="n"))
logAxis(2)

plot.oe.bet <- function(b,tm,...) {
	these.time <- tm[names(b),]
	b <- 1/b
	b <- log(b)
	# b <- 2*max(b) - b
	
	ext.complt <- these.time[,2] != min(these.time[,2])
	
	these.time <- these.time[ext.complt,]
	b <- b[ext.complt]
	
	b.cints <- classIntervals(b,10,intervalClosure="right")
	plot(b.cints,pal=tim.colors(10))
	
	b.col <- findColours(b.cints,tim.colors(10))
	plot(these.time[,1:2],col=b.col,xlim=rev(range(these.time[,1])),ylim=rev(range(these.time[,2])),...)
	abline(0,1)
}

par(mfcol=c(1,2))
plot.oe.bet(pbdb.sstat.ord.cor$beta,pbdb.ord.range)



##	relationship of beta with total richness
source("~/R_functions/logPlot.R")

pbdb.ord.max <- apply(pbdb.ord.div[,names(pbdb.sstat.ord.cor$beta)],2,max)

hist(pbdb.ord.max*pbdb.sstat.ord.cor$beta,breaks=12,prob=TRUE,
	 xlab=expression("(Maximum order diversity)*"*beta),ylab="Density",main="")
lines(density(pbdb.ord.max*pbdb.sstat.ord.cor$beta))

logPlot(pbdb.sstat.ord.cor$beta,pbdb.ord.max,log="xy",
		xlab=expression(beta),ylab="Maximum order diversity",
		mgp=c(2.5,1,0))

logPlot(pbdb.sstat.ord.cor$beta,pbdb.sstat.ord.cor$beta*pbdb.ord.max,log="xy",
		xlab=expression(beta),ylab=expression("(Maximum order diversity)*"*beta),
		mgp=c(2.5,1,0))

logPlot((pbdb.sstat.ord.cor$beta)^-1,pbdb.sstat.ord.cor$beta*pbdb.ord.max,log="xy",
		xlab=expression(sigma^2),ylab=expression("(Maximum order diversity)/"*sigma^2),
		mgp=c(2.5,1,0))

logPlot(pbdb.sstat.ord.cor$beta^-1,(pbdb.sstat.ord.cor$beta*pbdb.ord.max)^-1,log="xy",
		xlab=expression(sigma^2),ylab=expression(sigma^2*"/(Maximum order diversity)"),
		mgp=c(2.5,1,0))


plot.betRange(pbdb.sstat.ord.cor$beta*pbdb.ord.max^2,pbdb.ord.range,
			  xlab="Order lifetime",ylab=expression(beta*"*max.genera"^2),
			  yaxt="n",log="y",mgp=c(2.5,1,0))
logAxis(2)

plot.betTime(pbdb.sstat.ord.cor$beta*pbdb.ord.max^2,pbdb.ord.range,
			 frst=TRUE,mgp=c(2.5,1,0),
			 log="y",yaxt="n",
			 ylab=expression(beta*"*max.genera"^2))
logAxis(2)

