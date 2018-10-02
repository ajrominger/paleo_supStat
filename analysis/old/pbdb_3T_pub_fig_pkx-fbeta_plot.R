##  this version is FINAL, it makes p_k(x) and f(beta) plot
##  but still sources scripts and data in the 'Desktop/Research/paleo_supStat'
##  directory, that stuff should be migrated

##  to be migrated:
##      sstat_plotting_par.R
##      pbdb_3T_pub_lik_v2.R (and name change to 'pbdb_3T_pub_lik.R')
##      then need to migrate stuff that pbdb_3T_pub_lik_v2.R calls:
##           sstat_methods.R
##           make_pbdb_3Tpub.R (maybe need migrate from this too)

oldwd <- setwd('~/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub')

##	needed functions and data
source('../../code/sstat_plotting_par.R')
source("../../code/sstat_methods.R")

if(!exists("pbdb.sstat.ord.cor")) {
	cat('making sstat objects','\n')
	makePlot <- FALSE
	source("make_pbdb_3Tpub.R")
}

##  example taxa
eg.taxa <- c('Ammonitida', 'Proetida', 'Spiriferida')

##  extract data on example taxa
eg.taxa.traj <- lapply(eg.taxa, function(taxon) {
	this.traj <- pbdb.ord.div[,taxon]
	good <- (min(which(this.traj > 0)) - 1) : (max(which(this.traj > 0)) + 1)
	
	this.traj <- this.traj[good]
	this.time <- pbdb.time[rownames(pbdb.ord.div)][good]
	
	cbind(this.time, this.traj)
})


##  make eCDF for each order
all.pk.sub <- do.call(rbind,
                      lapply(pbdb.sstat.ord.cor$Px.sub, 
                             function(x) my.ecdf(scale(x), TRUE)))


##  begin plotting

with(plot.pars, {
	quartz(width=2.1*width, height=1.7*height)
	layout(matrix(c(1, 1, 2, 3), nrow=2, byrow=TRUE),
	       heights=c(0.75, 1))
	par(oma=mar, mgp=mgp)
	
	##  plot trajectories of example taxa
	par(mar=c(3, 0, 0, 0), xpd=NA)
	palette(c('firebrick', 'goldenrod', 'dodgerblue3'))
	paleoPlot(c(510, 20),
	          range(sapply(eg.taxa.traj, function(x) range(x[,2]))),
	          type='n', ylab='Standardized diversity')
	for(i in 1:length(eg.taxa)) {
		lines(eg.taxa.traj[[i]][,1] + 10*(i-1), eg.taxa.traj[[i]][,2], lwd=2, col=i)
	}
	
	##  plot distributions of all taxa
	par(mar=rep(0, 4))
	plot(all.pk.sub, col='gray',
	     xlab='Scaled fluctuations', ylab='Cumulative density')
	for(i in 1:length(eg.taxa)) {
		lines(my.ecdf(scale(pbdb.sstat.ord.cor$Px.sub[[eg.taxa[i]]]), TRUE), 
		      lwd=2, col=i)
	}
	curve(pnorm(x, lower.tail=FALSE), add=TRUE, lwd=2)
	
	##  plot f(beta)
	with(pbdb.sstat.ord.cor, {
		plot(my.ecdf(beta, TRUE), log='x', col='gray',
		     xaxt='n', yaxt='n', xlab=expression(beta), ylab='')
		curve(pgamma(x, gam.par[1], gam.par[2], lower.tail=FALSE),
		      lwd=2, add=TRUE)
	})
	logAxis(1)
})

dev.copy2pdf(file='~/Dropbox/Research/paleo_supStat/ms/science/version13/fig_pkx-fbeta.pdf')

setwd(oldwd)