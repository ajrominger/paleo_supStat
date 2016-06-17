##############  super-stat  ###############

##	control parameter for determining from where order-level data come
if(!exists("sepk2use")) sepk2use <- "sepk.cert"

##	control parameters for determining whether or not to compute different taxonomic scales
if(!exists("yesOrd")) yesOrd <- TRUE
if(!exists("yesCls")) yesCls <- TRUE
if(!exists("yesPhy")) yesPhy <- TRUE
if(!exists("yesOE")) yesOE <- FALSE

source("/Users/andrewrominger/Desktop/Research/paleo_supStat/sepk/sepk_obCreate.R")
sepk4ord <- get(sepk2use)	# either `sepkoski' or `sepk.cert'

##	needed packages
#library(MASS)
#
###	needed source code and data obs.
#load("sepk_sstat_obs.RData")
#source("~/R_functions/paleoPlot.R")
#source("~/R_functions/my_ecdf.R")
#source("~/R_functions/normLS.R")
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_comp.R")

##	get all diversity

##	create all the p_k(x|beta) for different metrics, taxonomic levels

# order-level
if(yesOrd) {
	sepk.ord <- split(sepk4ord[,5:ncol(sepk.cert)],sepk4ord$ord)
	sepk.ord.divflux <- make.flux.tab(sepk.ord)
	
	sepk.div <- apply(sapply(sepk.ord.divflux,function(x) x$div),1,sum)
	
	sepk.ord.dflux <- lapply(sepk.ord.divflux,function(x) event.flux(x$div))
	sepk.ord.oflux <- lapply(sepk.ord.divflux,function(x) event.flux(x$frst))
	sepk.ord.eflux <- lapply(sepk.ord.divflux,function(x) event.flux(x$last))

	sepk.ord.odflux <- lapply(sepk.ord.divflux,function(x) 
								{
									rate <- log(x$frst/sepk.div)
									rate <- rate[is.finite(rate)]
									event.flux(rate)
								})
	sepk.ord.edflux <- lapply(sepk.ord.divflux,function(x) {
									rate <- log(x$last/sepk.div)
									rate <- rate[is.finite(rate)]
									event.flux(rate)
								})
	
	##	diversity fluctuation superStat object
	sepk.sstat.ord.d <- sstat.comp(sepk.ord.dflux,minN=10,plotit=FALSE)
	
	##	origination/extinction fluctuation superStat objects
	if(yesOE) {
		sepk.sstat.ord.o <- sstat.comp(sepk.ord.oflux,minN=10,plotit=FALSE)
		sepk.sstat.ord.e <- sstat.comp(sepk.ord.eflux,minN=10,plotit=FALSE)
	}
}

# class-level
if(yesCls) {
	sepk.cls <- split(sepkoski[,5:ncol(sepkoski)],sepkoski$class)
	sepk.cls.divflux <- make.flux.tab(sepk.cls)
	
	sepk.cls.dflux <- lapply(sepk.cls.divflux,function(x) event.flux(x$div))
	
	##	diversity fluctuation superStat object
	sepk.sstat.cls.d <- sstat.comp(sepk.cls.dflux,minN=10,plotit=FALSE)

}

# phylum-level
if(yesPhy) {
	sepk.phy <- split(sepkoski[,5:ncol(sepkoski)],sepkoski$phy)
	sepk.phy.divflux <- make.flux.tab(sepk.phy)
	
	sepk.phy.dflux <- lapply(sepk.phy.divflux,function(x) event.flux(x$div))
	
	##	diversity fluctuation superStat object
	sepk.sstat.phy.d <- sstat.comp(sepk.phy.dflux,minN=10,plotit=FALSE)
}


#plot(my.ecdf(abs(unlist(sepk.sstat.ord.d$Px.sub)),TRUE),col="gray",log="xy")
#points(my.ecdf(abs(unlist(sepk.sstat.ord.d$Px.sub)),TRUE),type="s")
#curve(sepk.sstat.ord.d$PPx(x),col="red",add=TRUE)
#
#
#x <- sepk.sstat.ord.d
#
#plot(my.ecdf(unlist(x$Px.sub)))
#points(my.ecdf(unlist(y$Px.sub)),col=hsv(1,1,1,alpha=0.1))
#points(my.ecdf(unlist(sepk.sstat.ord.d$Px.sub)),col=hsv(0.5,1,1,alpha=0.2))
