############  order uniqueness for PBDB corrected  ##############

source("~/R_functions/ks_stat_pfun.R")
source("~/R_functions/den_fill.R")

##	depends on .../PBDB_3T_pub.R

##	function to make .flux from matrix
make.flux <- function(x,grp) {
	div <- aggregate(x,list(grp=grp),sum)
	
	apply(div[,-1],1,function(X)
					{
						flux <- diff(c(0,X))
						return(flux[flux != 0])
					})
}


##	matrix form makes easier to randomize, etc
pbdb.ma <- samp2site.spp(site=pbdb.dat$collections.10_my_bin,
						 spp=pbdb.dat$occurrences.genus_name,
						 abund=rep(1,nrow(pbdb.dat)))

pbdb.ma[pbdb.ma > 0] <- 1
pbdb.ma <- t(pbdb.ma)
pbdb.ma <- pbdb.ma[,order(pbdb.time,decreasing=TRUE)]

pbdb.ma.ord <- as.character(pbdb.dat$occurrences.order_name[match(rownames(pbdb.ma),pbdb.dat$occurrences.genus_name)])
pbdb.ma.cls <- as.character(pbdb.dat$occurrences.class_name[match(rownames(pbdb.ma),pbdb.dat$occurrences.genus_name)])

pbdb.ma.t3 <- pbdb.curv$Three.timer.sampling.stat[match(colnames(pbdb.ma),pbdb.curv$Bin.name)]
pbdb.ma.pub <- tapply(pbdb.dat$collections.reference_no,pbdb.dat$collections.10_my_bin,function(x) length(unique(x)))[colnames(pbdb.ma)]
pbdb.ma.pubEf <- pbdb.ma.pub^lm(log(ord.tbin.bias$T3.div)~log(ord.tbin.bias$tbin.pub))$coeff[2]

##	matrix of corrected occurrences (corrected by transitive property of 3T and resid)
pbdb.ma.cor <- t(t(pbdb.ma) * as.vector(pbdb.ma.t3 * pbdb.ma.pubEf)^-1)

##	checks out!!!
plot(scale(apply(pbdb.ma.cor,2,sum)),type="l",col="red")
lines(scale(rev(pbdb.samp$Mean.sampled.diversity)))

pbdb.ma.cor <- pbdb.ma.cor[,2:48]	# remove Cambrian 1 and Cenozoic 6

pbdb.flux.cor <- make.flux(pbdb.ma.cor,pbdb.ma.ord)
pbdb.sstat.ord <- sstat.comp(pbdb.flux.cor,minN=10)

plot.flux <- function(x,leg=TRUE,add.norm=TRUE,...) {
	pk <- x$raw.pk[x$incld]
	PPx <- x$PPx
	
	this.ecdf <- my.ecdf(abs(unlist(pk)),TRUE)
	
	logPlot(this.ecdf,log="xy",...)
	
	if(add.norm) {
		best.norm <- normLS(unlist(pk))$par
		curve(2*pnorm(x,0,best.norm[2],lower.tail=FALSE),col="blue",lwd=3,add=TRUE)
	}
	
#	curve(PPx(x),col="red",lwd=3,add=TRUE)
	
	if(leg) {
		legend("bottomleft",legend=c("Observed","Superstatistical","Best Gaussian"),
			   pch=1,lty=1,col=c("black","red","blue"),
			   lwd=c(0,1,1),pt.lwd=1,pt.cex=c(1,0,0),bty="n")
	}
}

quartz(width=4.5,height=4.5)
par(mar=c(5,4,3,1)+0.1)
plot.flux(pbdb.sstat.ord,xlab="|Fluctuations|",ylab="Cumulative density",leg=FALSE)
legend("bottomleft",legend=c("Observed data","Best normal"),col=c("black","blue"),
	   lty=1,pch=1,pt.cex=c(1,0),pt.lwd=c(1,0),lwd=c(0,2),bty="n")


##	permutation
nrun <- 100
pbdb.ord.permD <- numeric(nrun)

pbdb.ord.perm <- sample(pbdb.ma.ord)

##	begin permutation
for(i in 1:nrun) {
	##	permute order classification
	pbdb.ord.perm <- sample(pbdb.ord.perm)
	
	##	calc flux
	this.ord.flux <- make.flux(pbdb.ma.cor,pbdb.ord.perm)
	
	##	diversity fluctuation superStat object
	this.sstat <- sstat.comp(this.ord.flux,minN=10,plotit=FALSE)
	
	this.PPx <- function(x) 0.5 + 0.5*this.sstat$PPx(x,FALSE)
	
	pbdb.ord.permD[i] <- ks.stat.pfun(unlist(this.sstat$raw.pk),"this.PPx")
}

##	calc actual D
real.PPx <- function(x) 0.5 + 0.5*pbdb.sstat.ord$PPx(x,FALSE)
pbdb.ord.D <- ks.stat.pfun(unlist(pbdb.sstat.ord$raw.pk),"real.PPx")

quartz(width=4.5,height=4.5)
par(mar=c(5,4,3,1)+0.1)
den.fill(pbdb.ord.permD,xlim=c(0.05,0.09),xlab="D statistic",main="")

abline(v=pbdb.ord.D,col="red")
text(pbdb.ord.D+0.002,125,labels="Observed",pos=4,col="red")
arrows(pbdb.ord.D+0.003,125,pbdb.ord.D,110,col="red",length=0.1,angle=20)


##	calc sstat for classes
pbdb.sstat.cls <- sstat.comp(make.flux(pbdb.ma.cor,pbdb.ma.cls))
plot.sstat(pbdb.sstat.cls,leg=FALSE,add.norm=FALSE,xlab="|Fluctuations|",ylab="Cumulative density")
legend("bottomleft",legend=c("Observed","Superstatistical"),
	   pch=1,lwd=c(0,1),pt.lwd=c(1,0),pt.cex=c(1,0),col=c("black","red"),bty="n")