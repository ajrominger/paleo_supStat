#########  calculate beta and f(beta) for raw pbdb data  #########

##	load data
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb_read-in.R")

##	load functions
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_comp.R")
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_methods.R")


pbdb.ord.sstat <- sstat.comp(pbdb.ord.flux,minN=15,plotit=FALSE)


names(pbdb.ord.sstat)

##	object `sepk.sstat.ord.d' from .../sepk_supstat.R

sepk.col <- c(hsv(1,0.5,1),"red")
pbdb.col <- c(hsv(0.6,0.5),"blue")


##	f(beta)
logPlot(my.ecdf(sepk.sstat.ord.d$beta,TRUE),xlim=range(sepk.sstat.ord.d$beta,pbdb.ord.sstat$beta),
		log="x",col=sepk.col[1],xlab=expression(beta),ylab="Cummulative density",cex.lab=1.2)
points(my.ecdf(pbdb.ord.sstat$beta,TRUE),col=pbdb.col[1])

curve(with(sepk.sstat.ord.d,pgamma(x,gam.par[1],gam.par[2],lower.tail=FALSE)),col=sepk.col[2],add=TRUE)
curve(with(pbdb.ord.sstat,pgamma(x,gam.par[1],gam.par[2],lower.tail=FALSE)),col=pbdb.col[2],add=TRUE)

##	with CI
pbdb.beta.cdf <- my.ecdf(pbdb.ord.sstat$beta,TRUE)
pbdb.var.n <- sapply(pbdb.ord.sstat$Px.sub[pmatch(pbdb.beta.cdf[,1],pbdb.ord.sstat$beta)],length)
pbdb.beta.ci <- 1/cbind((1/pbdb.beta.cdf[,1])*qchisq(0.025,pbdb.var.n-1)/(pbdb.var.n-1),
						(1/pbdb.beta.cdf[,1])*qchisq(0.975,pbdb.var.n-1)/(pbdb.var.n-1))

logPlot(pbdb.beta.cdf,xlim=c(0.0001,1),
		log="x",col=pbdb.col[2],xlab=expression(beta),ylab="Cummulative density",cex.lab=1.2)


segments(x0=pbdb.beta.ci[,1],x1=pbdb.beta.ci[,2],y0=pbdb.beta.cdf[,2],col=pbdb.col[1])

##	P(x)
logPlot(my.ecdf(abs(unlist(sepk.sstat.ord.d$Px.sub)),TRUE),
		log="xy",col=sepk.col[1],xlab="|Fluctuations|",ylab="Cummulative density")
points(my.ecdf(abs(unlist(pbdb.ord.sstat$Px.sub)),TRUE),col=pbdb.col[1])

curve(sepk.sstat.ord.d$PPx(x),col=sepk.col[2],add=TRUE)
curve(pbdb.ord.sstat$PPx(x),col=pbdb.col[2],add=TRUE)

logPlot(my.ecdf(abs(unlist(pbdb.ord.sstat$Px.sub)),TRUE),log="xy",col="red",
		xlab="|Absolute fluctuations|",ylab="Cummulative density",main="P(x) for PBDB")
curve(pbdb.ord.sstat$PPx(x),add=TRUE)
legend("bottomleft",legend=c("Observed data","Super-statistic theory"),
	   col=c("red","black"),pch=c(1,NA),lty=c(NA,1),bty="n")

pbdb.beta.share <- pbdb.ord.sstat$beta[toupper(names(pbdb.ord.sstat$beta)) %in% names(sepk.sstat.ord.d$beta)]
sepk.beta.share <- sepk.sstat.ord.d$beta[names(sepk.sstat.ord.d$beta) %in% toupper(names(pbdb.ord.sstat$beta))]

all(names(sepk.beta.share) == toupper(names(pbdb.beta.share)))  # make sure all in same order

plot(sepk.beta.share,pbdb.beta.share,
	 xlab=expression(beta~"from Sepkoski"),
	 ylab=expression(beta~"from PBDB"))

pbdb.beta.out <- matrix(cbind(sepk.beta.share,pbdb.beta.share)[pbdb.beta.share > 0.4,],ncol=2)
rownames(pbdb.beta.out) <- names(pbdb.beta.share)[pbdb.beta.share > 0.4]

sepk.beta.out <- cbind(sepk.beta.share,pbdb.beta.share)[sepk.beta.share > 0.4,]
rownames(sepk.beta.out) <- names(pbdb.beta.share)[sepk.beta.share > 0.4]

text(pbdb.beta.out,labels=rownames(pbdb.beta.out),pos=1)
text(sepk.beta.out,labels=rownames(sepk.beta.out),pos=2)



plot(sepk.beta.share,pbdb.beta.share,xlim=c(0,0.3),ylim=c(0,0.3),
	 xlab=expression(beta~"from Sepkoski"),
	 ylab=expression(beta~"from PBDB"))

