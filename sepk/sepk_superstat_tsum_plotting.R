###########  sepk tsum plotting  ###########

##	load simulated data--!!!!!!!!!need improve!!!!!!!!!, missing `sepk.ord.divflux' etc....
load("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_sstat_tsum_obs.RData")

##	function to summarize simulation output.
##	specifically, interpolate values in the eCDF. assumes log x-axis
ecdf.intrp <- function(x,min1=TRUE,n=100,alpha=0.05,ymin=0.001) {
	x.rng <- range(sapply(x,function(X) range(X[,1])))
	if(min1) x.rng[1] <- 1
	
	x.new <- exp(seq(from=log(x.rng[1]),to=log(x.rng[2]),length=n-1))
	x.new <- c(x.new,exp(log(x.new[n-1]) + diff(log(x.new[1:2]))))
	
	y.appx <- sapply(x,function(X) approx(X,xout=x.new,method="constant",yright=ymin,rule=2,f=1)$y)
	
#	y.ci <- apply(y.appx,1,quantile,probs=c(alpha/2, 0.5, 1 - alpha/2))
	y.ci <- apply(y.appx,1,function(x) c(quantile(x,alpha/2),mean(x),quantile(x,1 - alpha/2)))
	
	return(cbind(x=x.new, ylow=y.ci[1,], ymid=y.ci[2,], yhi=y.ci[3,]))
}

##	funciton to plot eCDF summary
plot.ecdf.intrp <- function(X,dat,...) {
	logPlot(dat,xlim=range(X[,1]),ylim=range(X[,-1],dat[,2],na.rm=TRUE),...)
	for(i in 2:4) points(X[,c(1,i)],lty=ifelse(i %% 2 == 0, 2, 1),type="l",col="red")
}

##	function to compute tail exponent
tsum.tail <- function(x,tail.n=20) {
	x <- log(x[1:tail.n,])
	expo <- lm(x[,2]~x[,1])$coeff[2]
	
	return(as.numeric(expo))
}


##	CDF
sepk.sstat.ordQ <- apply(sapply(sepk.ord.divflux,function(x) x$flux[-recent.rm]),1,sum)
Qx.real <- my.ecdf(abs(scale(sepk.sstat.ordQ,center=FALSE,scale=FALSE)),TRUE)

Qx.sim.intrp <- ecdf.intrp(tsum.sim,ymin=0.001)

Qx.sim.intrp[Qx.sim.intrp < min(Qx.real[,2])] <- NA
intrp.stop <- min(which(apply(Qx.sim.intrp[,-1],1,function(x) all(is.na(x)))))-1
if(intrp.stop < nrow(Qx.sim.intrp)) Qx.sim.intrp <- Qx.sim.intrp[1:intrp.stop,]

plot.ecdf.intrp(Qx.sim.intrp,Qx.real,log="xy",
				xlab="|Global fluctuations|",ylab="Cummulative density")

legend("bottomleft",legend=c("Observed","Superstatistics"),
	   col=c("black","red"),pch=c(1,22),pt.bg=hsv(alpha=0.3),pt.cex=1:2,
	   lty=1,lwd=0:1,pt.lwd=1:0,bty="n")


#sepk.ord.fluxXX <- sepk.ord.divflux[[1]][1:73,c(1,3)][order(abs(sepk.sstat.ordQ),decreasing=TRUE),][1:5,]
#sepk.ord.fluxXX[,2] <- sepk.sstat.ordQ[order(abs(sepk.sstat.ordQ),decreasing=TRUE)][1:5]
#sepk.ord.fluxXX
#
#
###	tail exponent
#Qx.sim.tail <- sapply(tsum.sim[1:192],tsum.tail)
#Qx.act.tail <- tsum.tail(my.ecdf(abs(scale(sepk.sstat.ordQ,center=FALSE,scale=FALSE)),TRUE))
#
#den.fill(Qx.sim.tail,xlab="Tail exponent",main="")
#abline(v=Qx.act.tail,col="red")
#
#
###	sigma of each order
#tsum.beta.sum <- t(apply(tsum.sim.beta,1,quantile,probs=c(0.025,0.5,0.975)))
#plot(sepk.time.orig.sig[order(sepk.time.orig.sig[,3]),3]^-2,tsum.beta.sum[order(sepk.time.orig.sig[,3]),3],type="l")
#points(sepk.time.orig.sig[order(sepk.time.orig.sig[,3]),3]^-2,tsum.beta.sum[order(sepk.time.orig.sig[,3]),1],type="l")
#points(sepk.time.orig.sig[order(sepk.time.orig.sig[,3]),3]^-2,tsum.beta.sum[order(sepk.time.orig.sig[,3]),2],type="l")
#abline(0,1,lty=2)
#
#
#
#fbeta.sim <- ecdf.intrp(apply(tsum.sim.beta,2,my.ecdf,complement=TRUE),min1=FALSE,ymin=0)
#
#plot.ecdf.intrp(fbeta.sim,my.ecdf(sepk.time.orig.sig[,3]^-2,TRUE),log="x",col="red",cex=0.5,pch=16)
#
