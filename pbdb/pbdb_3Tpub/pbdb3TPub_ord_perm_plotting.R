###############  plot results from ord-perm for 3TPub corrected PBDB  ###############

oldwd <- setwd('~/Dropbox/research/paleo_supStat/pbdb/pbdb_3Tpub')

##	function for plotting desnities
source("~/R_functions/den_fill.R")
source("~/R_functions/ks_stat_pfun.R")

##	depends on .../make_pbdb_3Tpub.R
if(!(exists("pbdb.ord.div") & exists("pbdb.gen.occ3TP") & exists("pbdb.sstat.ord.cor"))) {
	makePlot <- FALSE
	source("make_pbdb_3Tpub.R")
}

##	load permuted d-stat
load("ordperm_dstat.RData")

pbdb.ord.dstat <- with(pbdb.sstat.ord.cor,ks.stat.pfun(unlist(raw.pk),function(x) 0.5 + 0.5*PPx(x,FALSE)))
pbdb.cls.dstat <- with(pbdb.sstat.cls.cor,ks.stat.pfun(unlist(raw.pk),function(x) 0.5 + 0.5*PPx(x,FALSE)))
pbdb.phy.dstat <- NA	# STUB!!!

##	plot it!
source('../../code/sstat_plotting_par.R')

with(plot.pars, {
	quartz(width=width,height=height)
	par(mar=mar,mgp=mgp)
})

den.fill(ordperm.dstat,alpha=0.05,
		 xlim={
		 	x <- range(ordperm.dstat,pbdb.ord.dstat,pbdb.cls.dstat,pbdb.phy.dstat,na.rm=TRUE)
		 	x + c(-1,1)*0.04*max(x)
		 },
		 col="gray20",xlab="D-statistic",main="")
abline(v=c(pbdb.ord.dstat,pbdb.cls.dstat,pbdb.phy.dstat),
	   # col=hsv(s=c(1,0.4,0.4)),
	   lwd=3,lty=1:2)

text(c(pbdb.ord.dstat,pbdb.cls.dstat,pbdb.phy.dstat),rep(par("usr")[4]*0.8,3),
	 labels=c("Orders","Classes","Phyla"),srt=90,
	 # col=hsv(s=c(1,0.3,0.4)),
	 adj=c(0.5,1.5))

# mtext(date(),side=4,line=0.1)
# mtext("``pbdb3TP_ord_perm_plotting.R''",side=4,line=1.2)

setwd(oldwd)