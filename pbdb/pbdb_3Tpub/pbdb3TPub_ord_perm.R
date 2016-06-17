###############  ord-perm for 3TPub corrected PBDB  ###############

oldwd <- setwd("~/Desktop/Research/paleo_supStat/pbdb/pbdb_3Tpub")

##	depends on .../make_pbdb_3Tpub.R
if(!(exists("pbdb.ord.div") & exists("pbdb.gen.occ3TP") & exists("pbdb.sstat.ord.cor"))) {
	makePlot <- FALSE
	source("make_pbdb_3Tpub.R")
}

##	functions for d-stat
ks.stat.pfun <- function(x,pfun,...) {
	n <- length(x)
	y <- pfun
	
	x <- y(sort(x),...) - (0:(n - 1))/n
	stat <- max(c(x, 1/n - x))
	
	return(stat)
}

##	convenience function to split up time-bin by genus occurance matrix
##	into order flux list
genOcc2ord <- function(genOcc,genMap,ords) {
	lapply(ords,function(x) {
		raw.div <- rowSums(genOcc[,genMap==x,drop=FALSE])
		raw.flux <- diff(c(0,raw.div))
		return(raw.flux[raw.flux != 0])
	})
}

##	vector mapping genera (columns) in pbdb.gen.occ3TP to orders
gen2ord <- as.character(with(pbdb.dat,occurrences.order_name[match(colnames(pbdb.gen.occ3TP),occurrences.genus_name)]))
all.ord <- colnames(pbdb.ord.div)

##	begin permutation
nsim <- 1000
ordperm.dstat <- replicate(nsim,{
	new.gen2ord <- sample(gen2ord)
	new.flux <- genOcc2ord(pbdb.gen.occ3TP,new.gen2ord,all.ord)
	
	this.sstat <- sstat.comp(new.flux,minN=15,plotit=FALSE)
	
	this.PPx <- function(x) 0.5 + 0.5*this.sstat$PPx(x,FALSE)
	ks.stat.pfun(unlist(this.sstat$raw.pk),this.PPx)
	cat(proc.time()[1],'\n')
})

##	save it
save(ordperm.dstat,file="ordperm_dstat.RData")

setwd(oldwd)