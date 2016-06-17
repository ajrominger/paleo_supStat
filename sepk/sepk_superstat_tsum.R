##	depends on objects made in .../sepk_supstat.R
yesCls <- FALSE
yesPhy <- FALSE
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_supstat.R")

source("~/R_functions/my_ecdf.R")
source("~/R_functions/logPlot.R")
source("~/R_functions/den_fill.R")

##	function to generate p_k(x|beta). output is a matrix,
##	rows = time, cols = time series of fluxes
pk.gen <- function(max.time,		# max duration of time
				   time.orig.sig) {	# matrix of durration, origination and standard deviation for each taxon
	
	out <- apply(time.orig.sig,1,function(x)
					{
						fseries <- csum.pos.norm(x[1],x[3],ext=(sum(x[1:2]) - 1 < max.time))
						temp <- numeric(max.time+1)  # last element will be beta
						temp[(1:x[1]) + x[2] - 1] <- fseries$series
						temp[max.time+1] <- fseries$beta
						return(temp)
						
					})
	
	return(list(beta=out[max.time+1,],series=out[-(max.time+1),]))
}

##	function generates one series of normal fluctuations whoes cumulative sum
##	is alway greater than equal to 0
csum.pos.norm <- function(n,sd,ext,ntry=1000) {
	this.try <- csum.pos.test(n,sd,ext,ntry)
	
	if(nrow(this.try) == 0) {
		for(i in 1:ntry) {
			this.try <- csum.pos.test(n,sd,ext,ntry)
#			print(i)
			if(nrow(this.try) != 0) break;
		}
	}
	
	out <- this.try[,sample(ncol(this.try),1)]
	
	return(list(beta=1/var(out),series=out))
}

##	function generates an attempt to make several normal series whoes cumulative sum
##	is alway greater than equal to 0. output is a matrix, columns are the series
csum.pos.test <- function(n,sd,ext,ntry) {
	X <- replicate(ntry*2,{
								x <- round(rnorm(n,sd=sd))
								end.cond <- ifelse(ext, sum(x) == 0, sum(x) > 0)
								if(all(cumsum(x[-n]) >= 0) & end.cond) {
									return(x)
								} else {
									return(numeric(0))
								}
							},simplify=FALSE)
	X <- do.call(cbind,X)
	
	return(X)
}

csum.pos.norm(50,0.30,100)
range(sepk.sstat.ord.d$beta^-0.5)


###############  simulations  ###############

##	compute time.orig.sig from sepk.ord.divflux and sepk.sstat.ord.d
##	extract gamma pars, ntaxa,...
##	in all cases, remove Recent = time 74
recent.rm <- 74
sepk.time.orig.sig <- t(sapply(sepk.ord.divflux, #[sepk.sstat.ord.d$incld],
							 function(x)
							 	{
							 		these.occ <- which(x$div[-recent.rm] > 0)
							 		return(c(diff(range(these.occ))+1,min(these.occ)))
							 	}))
sepk.time.orig.sig <- cbind(sepk.time.orig.sig,NA)
# those only in Recent will give -Inf
sepk.time.orig.sig <- sepk.time.orig.sig[is.finite(sepk.time.orig.sig[,1]),]

# take sd from fitted vals
#sepk.time.orig.sig[,3] <- sepk.sstat.ord.d$beta^-0.5
# make sure match up
#all(names(sepk.sstat.ord.d$beta^-0.5) == rownames(sepk.time.orig.sig))

nrun <- 500
ntaxa <- nrow(sepk.time.orig.sig)
max.time <- nrow(sepk.ord.divflux[[1]][-recent.rm,])
gamma.shape <- sepk.sstat.ord.d$gam.par[1]
gamma.rate <- sepk.sstat.ord.d$gam.par[2]


tsum.sim.raw <- tsum.sim <- vector("list",nrun)
tsum.sim.beta <- matrix(NA,ntaxa,nrun)


for(i in 1:nrun) {
	# generate standard devs
	sepk.time.orig.sig[,3] <- rgamma(ntaxa,gamma.shape,gamma.rate)^-0.5
	
	these.pk <- try(pk.gen(max.time,sepk.time.orig.sig))
	
	if(class(these.pk) != "try-error") {
		# one series of global flux
		this.series <- apply(these.pk$series,1,sum)
		this.series <- this.series[this.series != 0]   # ignore 0 just like in original analysis of P(x)
	} else {
		this.series <- NA
		these.pk <- list(beta=NA)
	}
#	print(all(cumsum(this.series) >= 0))
	
	# raw simulated data
	tsum.sim.raw[[i]] <- this.series
	
	# cdf
	tsum.sim[[i]] <- my.ecdf(abs(this.series),complement=TRUE)
	tsum.sim.beta[,i] <- these.pk$beta
	print(i)
}
tsum.sim <- tsum.sim[1:21]
##	get rid of error simulations
tsum.sim <- tsum.sim[sapply(tsum.sim,nrow) > 0]

save(tsum.sim.raw,tsum.sim,sepk.ord.divflux,sepk.sstat.ord.d,recent.rm,
	 file="/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_sstat_tsum_obs.RData")
