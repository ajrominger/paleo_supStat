##  data and sstat analysis
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub/make_pbdb_3Tpub.R")

##	function to generate a normal flux trajectory that originates and goes extint at the appropriate time

pk.gen <- function(frst,last,tot.rng,sig,ntry=1000) {
	out <- numeric(diff(tot.rng)+1)  # vector to be filled with fluctuations
	
	ext <- last < tot.rng[2]  # is the lineage extint?
	
	n <- last - frst + 1  # number of draws from normal
	
	this.gaus <- cand.gaus(n,ext,sig,ntry)  # a guassian flux trajectory
	
	i <- 1
	while(i < 4 & any(is.na(this.gaus))) {  # in case trajectory is unacceptable
		print(i)
		this.gaus <- cand.gaus(n,ext,sig,ntry)
	}
	
	out[frst:last] <- this.gaus
	return(out)
}

##  produces canditate trajectories for pk.gen
cand.gaus <- function(n,ext,sig,ntry) {
	these.gaus <- replicate(ntry,cumsum(rnorm(n,sd=sig)))  # candidate trajectories
	
	sign.these.gaus <- these.gaus > 0
	
	if(ext) {  # if lineage is extinct sum at n must be <= 0
		these.ok <- which(apply(sign.these.gaus[-n,,drop=FALSE],2,all) & !sign.these.gaus[n,])
	} else {   # if not extinct then must be > 0
		these.ok <- which(apply(sign.these.gaus,2,all))
	}
	
	if(length(these.ok) < 1) {
		return(NA)
	} else {
		return(diff(c(0,these.gaus[,sample(these.ok,1)])))
	}
}


##  calculate first/lasts from PBDB
pbdb.frst.last <- t(apply(pbdb.ord.div,2,function(x) range(which(x > 0)) + c(0,1)))
pbdb.rng <- c(1,nrow(pbdb.ord.div)+1)

##  simulate Q(x)
nrun <- 500
sim.Q.raw <- matrix(NA,nrow=pbdb.rng[2],ncol=nrun)

for(i in 1:nrun) {
	print(paste("*~*~*~*~*~*ITTERATION:",i))
	this.sim <- apply(pbdb.frst.last,1,function(x)
						pk.gen(x[1],x[2],pbdb.rng,
							   1/sqrt(rgamma(1,pbdb.sstat.ord.cor$gam.par[1],pbdb.sstat.ord.cor$gam.par[2])))
						)
	sim.Q.raw[,i] <- apply(this.sim,1,sum)
}

##  plotting
load("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub/Qsim.RData")
load("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub/pbdb_ordDiv_3Tpub.RData")

##  convenience function for interpolating
ecdf.interp <- function(x,knots) {
	this.ecdf <- my.ecdf(abs(x),complement=TRUE)
	
	approx(this.ecdf,xout=knots,method="constant",yleft=1,yright=0,f=1)$y
}

##  values to interpolate over
sim.Q.xval <- exp(seq(log(min(abs(sim.Q.raw))),log(max(abs(sim.Q.raw))),by=0.1))

##  interpolated CI for CDF
sim.Q.interp <- apply(sim.Q.raw,2,ecdf.interp,knots=sim.Q.xval)
sim.Q.ci <- apply(sim.Q.interp,1,function(x) c(mean=mean(x),quantile(x,prob=c(0.025,0.975))))

##  observed CDF
pbdb.Q.cor <- my.ecdf(abs(diff(c(0,apply(pbdb.ord.div,1,sum)))),TRUE)

##  clip sim.Q.ci to range of observed
this.ylim <- c(min(pbdb.Q.cor[,2])*0.9,1)
sim.Q.ci[sim.Q.ci < this.ylim[1]] <- this.ylim[1]
sim.Q.ci <- sim.Q.ci[,!apply(sim.Q.ci,2,function(x) sum(duplicated(x))>1)]
sim.Q.xval2 <- sim.Q.xval[1:ncol(sim.Q.ci)]

quartz(width=4,height=4)
par(mar=c(4,4,2,1)+0.1,mgp=c(2.5,1,0))
logPlot(pbdb.Q.cor,xlim=range(sim.Q.xval2),ylim=this.ylim,
		xlab="|Fluctuations|",ylab="Cumulative density",main="Q(x) for PBDB",
		log="xy",type="n")

polygon(x=c(sim.Q.xval2,rev(sim.Q.xval2)),
		y=c(sim.Q.ci[2,],rev(sim.Q.ci[3,])),
		col="gray85",border=NA)
points(my.ecdf(abs(diff(c(0,apply(pbdb.ord.div,1,sum)))),TRUE))
lines(sim.Q.xval2,sim.Q.ci[1,],col="gray55",lwd=1.5)

legend("bottomleft",legend=c("Observed data","Simulation mean","Simulation 95% CI"),
	   pch=c(1,1,15),lwd=c(0,1.5,0),pt.cex=c(1,0,1.5),pt.lwd=c(1,0,1),
	   col=c("black","gray55","gray85"),bty="n")
