#############  makes figures for sepkoski supplament  #############

##	needed functions
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_methods.R")
source("~/R_functions/den_fill.R")

##	sstat objects (ord, cls, phy)
yesOrd <- TRUE
yesCls <- TRUE
yesPhy <- TRUE
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_supstat.R")

##	D-statistic permutation
load("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_sstat_ord_perm.RData")

##	directory to save figs
save.dir <- "/Users/andrewrominger/Desktop/Research/paleo_supStat/ms/current_version"

##	P(x) fig for orders, classes and phyla
quartz(width=8,height=8/3)
layout(matrix(1:3,nrow=1))

par(mar=c(3.5,3,1,0)+0.1,mgp=c(2,0.75,0),cex.lab=1.3)
plot(sepk.sstat.ord.d,xlab="",ylab="Cumulative density")

par(mar=c(3.5,1.8,1,1.2)+0.1,mgp=c(2.5,1,0))
plot(sepk.sstat.cls.d,xlab="|Global fluctuations|",ylab="",add.legend=FALSE)

par(mar=c(3.5,0.6,1,2.4)+0.1,mgp=c(2.5,1,0))
plot(sepk.sstat.phy.d,xlab="",ylab="",add.legend=FALSE)

dev.print(pdf,file=paste(save.dir,"figSupp_sepk_P(x).pdf",sep="/"),width=8,height=8/3)


##	D-statistic
quartz(width=4,height=4)
par(mar=c(3.5,3.5,1,0.5)+0.1,mgp=c(2.5,1,0))
den.fill(sepk.ord.permD,xlim=range(0.8*sepk.actl.D,1.05*sepk.ord.permD),
		 xlab="D-statistic",main="",frame.plot=FALSE)
abline(v=sepk.actl.D,col="red",lwd=2)
legend("topleft",legend=c("Permuted orders","Observed"),
	   lwd=c(0,2),pch=15,pt.cex=c(2,0),col=c("gray","red"),
	   bg="white",box.col="transparent")
box()

dev.print(pdf,file=paste(save.dir,"figSupp_sepk_d-stat.pdf",sep="/"),width=4,height=4)


##	Q(x) sim
quartz(width=4,height=4)
par(mar=c(3.5,3.5,1,0.5)+0.1,mgp=c(2.5,1,0))
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_superstat_tsum_plotting.R")

dev.print(pdf,file=paste(save.dir,"figSupp_sepk_Qx_sim.pdf",sep="/"),width=4,height=4)