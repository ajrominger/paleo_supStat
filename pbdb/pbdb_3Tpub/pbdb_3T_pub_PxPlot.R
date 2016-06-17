oldwd <- setwd('~/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub')
source('../../code/sstat_plotting_par.R')

with(plot.pars, {
	quartz(width=width*2.1,height=height)
	par(oma=mar, mgp=mgp, mfrow=c(1,2), mar=c(0,0,0,0.5))
	# browser()
	source('pbdb_3T_pub_lik_v2.R')
	text(10^par("usr")[2], 10^par("usr")[4], labels="A", adj=adj, cex=cex.txt)
	ylim <- 10^par('usr')[3:4]
	
	par(mar=c(0,0.5,0,0))
	plot(pbdb.sstat.cls.cor,add.legend=FALSE,y.lim=ylim,yaxs='i',yaxt='n')
	text(10^par("usr")[2], 10^par("usr")[4], labels="B", adj=adj, cex=cex.txt)
	
	mtext("Cumulative density",side=2,line=2,outer=TRUE)
	mtext("|Fluctuations|",side=1,line=2,outer=TRUE)
})

setwd(oldwd)