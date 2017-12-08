##  this version is FINAL, it makes P(x) plot (both for main text and supp)
##  but still sources scripts and data in the 'Desktop/Research/paleo_supStat'
##  directory, that stuff should be migrated

##  to be migrated:
##      sstat_plotting_par.R
##      pbdb_3T_pub_lik_v2.R (and name change to 'pbdb_3T_pub_lik.R')
##      then need to migrate stuff that pbdb_3T_pub_lik_v2.R calls:
##           sstat_methods.R
##           make_pbdb_3Tpub.R (maybe need migrate from this too)

oldwd <- setwd('~/Dropbox/Research/paleo_supStat/pbdb/pbdb_3Tpub')
source('../../code/sstat_plotting_par.R')

##  plot P(x) for orders with likelihood CI's
with(plot.pars, {
	# quartz(width=width,height=height)
	par(mar=mar, mgp=mgp)
	# browser()
    source('pbdb/pbdb_3Tpub/pbdb_3T_pub_lik_v2.R')
	
	mtext("Cumulative density",side=2,line=2,outer=TRUE)
	mtext("|Fluctuations|",side=1,line=2,outer=TRUE)
})

dev.copy2pdf(file='~/Dropbox/Research/paleo_supStat/ms/science/version13/fig_Px.pdf')

##  plot P(x) for classes
with(plot.pars, {
	quartz(width=width,height=height)
	par(mar=mar, mgp=mgp)
	# browser()
	source('pbdb_3T_pub_lik_v2.R')
	
	plot(pbdb.sstat.cls.cor,add.legend=TRUE, xlab='', ylab='')
	
	mtext("Cumulative density",side=2,line=2)
	mtext("|Fluctuations|",side=1,line=2)
})

dev.copy2pdf(file='~/Dropbox/Research/paleo_supStat/ms/science/supp/figSupp_Px_cls.pdf')


setwd(oldwd)