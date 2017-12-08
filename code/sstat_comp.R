source('~/R_functions/my_ecdf.R')
source("~/R_functions/normLS.R")
source("~/R_functions/gammaLS.R")
source("~/R_functions/logPlot.R")
source('~/Dropbox/Research/paleo_supStat/code/Px_gam.R')

library(distr)

sstat.comp <- function(grp.data,minN=15,xlab="Absolute Fluctuation",ylab="Cumulative Density",leg=TRUE,plotit=TRUE) {
	these2use <- sapply(grp.data,length) >= minN
	p2use <- grp.data[these2use]
	
	cat("computing Gaussian fit for p_k(x|sigma) \n")
	pk.par <- sapply(p2use,function(x) unlist(normLS(x)[c("par","value")]))
	pk.par <- t(pk.par)
	colnames(pk.par) <- c("mu","sig","ss")
	
	cat("re-centering \n")
	for(i in 1:length(p2use)) {
		p2use[[i]] <- p2use[[i]] - pk.par[i,"mu"]
	}
	
	cat("computing f(beta) \n")
	f.beta.par <- gammaLS(1/(pk.par[,"sig"])^2)$par
	fuent.par <- c(n=2*f.beta.par[1],b0=f.beta.par[1]*f.beta.par[2])
	
	cat("computing P(x) \n")
	this.Px <- function(x) Px.gam(x,f.beta.par[1],f.beta.par[2])
	this.PPx <- function(x,comp=TRUE) PPx.gam(x,f.beta.par[1],f.beta.par[2],comp)
	
	##	functions using library `distr'
	sstat.dist <- "stub"#AbscontDistribution(d=this.Px)
	dsstat <- "stub"#function(x) sstat.dist@d (x)
	psstat <- "stub"#function(q,lower.tail=FALSE,log.p=FALSE) sstat.dist@p(q,lower.tail=lower.tail,log.p=log.p)
	qsstat <- "stub"#function(p,lower.tail=TRUE,log.p=FALSE) sstat.dist@q(p,lower.tail=lower.tail,log.p=log.p)
	rsstat <- "stub"#function(n) sstat.dist@r(n)
	
	if(plotit) {
		cat("plotting \n")
		this.ecdf <- my.ecdf(abs(unlist(p2use)),TRUE)
		
		logPlot(this.ecdf,xlab=xlab,ylab=ylab,log="xy",type="p")
		points(my.ecdf(abs(unlist(grp.data)),TRUE),col=hsv(0,1,0,alpha=0.4))
		
		curve(this.PPx(x),col="red",add=TRUE)
		
		# best.norm <- normLS(unlist(p2use))$par
		best.norm <- c(mean(unlist(p2use)),sd(unlist(p2use)))
		
		curve(2*pnorm(x,best.norm[1],best.norm[2],lower.tail=FALSE),col="blue",add=TRUE)
		
		if(leg) {
			legend("bottomleft",legend=c("All data",paste("Data w/ N >=",minN),"Superstatistical","Best Gaussian"),
				   pch=1,lty=1,col=c(hsv(0,1,0,alpha=0.4),"black","red","blue"),
				   lwd=c(0,0,1,1),pt.lwd=1,pt.cex=c(1,1,0,0),bty="n")
		}
	}
	
	out <- list(gam.par=f.beta.par,sspar=fuent.par,beta=1/(pk.par[,"sig"])^2,sumSq=pk.par[,"ss"],minN=minN,raw.pk=p2use,
				Px.raw=grp.data,Px.sub=p2use,incld=these2use,Px=this.Px,PPx=this.PPx,
				dsstat=dsstat,psstat=psstat,qsstat=qsstat,rsstat=rsstat)
	
	class(out) <- "sstat"
	
	return(out)
}

