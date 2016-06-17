#############  permatuation test of order-level uniqueness  #############

##	load data and sstat obs
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_supstat.R")
sepk2use <- "sepk.cert"
sepk4ord <- get(sepk2use)	# either `sepkoski' or `sepk.cert'

##	needed source code and data obs.
source("~/R_functions/ks_stat_pfun.R")
#source("~/R_functions/den_fill.R")

##	each iteration, the order tag will be permuted
sepk.ord.perm <- sepk4ord
#sepk.ord.perm$ord <- as.character(sepk.ord.perm$ord)

##	simulation parameters
nrun <- 499
sepk.ord.permD <- numeric(nrun)

##	begin permutation
for(i in 1:nrun) {
	##	permute order classification
	sepk.ord.perm$ord <- sample(sepk.ord.perm$ord)
	
	##	calc flux
	this.sepk.ord <- split(sepk.ord.perm[,5:ncol(sepk.cert)],sepk.ord.perm$ord)
	this.ord.divflux <- make.flux.tab(this.sepk.ord)
	
	this.ord.dflux <- lapply(this.ord.divflux,function(x) event.flux(x$div))
	
	##	diversity fluctuation superStat object
	this.sstat <- sstat.comp(this.ord.dflux,minN=10,plotit=FALSE)
	
	this.PPx <- function(x) 0.5 + 0.5*this.sstat$PPx(x,FALSE)
	
	sepk.ord.permD[i] <- ks.stat.pfun(unlist(this.sstat$raw.pk),"this.PPx")
}

##	real D-statistic
sepk.actl.PPx <- function(x) 0.5 + 0.5*sepk.sstat.ord.d$PPx(x,FALSE)
sepk.actl.D <- ks.stat.pfun(unlist(sepk.sstat.ord.d$raw.pk),"sepk.actl.PPx")

#den.fill(sepk.ord.permD,xlim=c(0.15,0.33));abline(v=sepk.actl.D,col="red")

save(sepk.actl.D,sepk.ord.permD,file="/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_sstat_ord_perm.RData")