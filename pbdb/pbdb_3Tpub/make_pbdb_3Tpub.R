#######  corrects diversity by 3-timer then number of pubs  #######
if(!exists("makePlot")) makePlot <- FALSE


oldwd <- setwd('~/Dropbox/Research/paleo_supStat')

##  convenience function to produce a matrix of time by ord with cells
##  of corrected diversity
source('code/pbdb_3t_pub.R')

##	load other needed funcitons
paleoPlot <- function(...) {
    plot(...)
    socorro::paleoAxis(1)
}
# source("~/R_functions/paleoPlot.R")

samp2site.spp <- socorro::tidy2mat
# source("~/R_functions/samp2site_spp.R")

logPlot <- function(..., log) {
    plot(..., log = log)
    s <- (1:2)[c(grepl('x', log), grepl('y', log))]
    
    socorro::logAxis(s)
}
# source("~/R_functions/logPlot.R")

my.ecdf <- socorro::simpECDF
# source("~/R_functions/my_ecdf.R")

source('code/sstat_comp.R')
source('code/sstat_methods.R')

##########  load data  ##########
setwd('data/old/pbdb_2013-05-28')

##	raw occurence data
pbdb.dat <- read.csv("marInv-occs.csv")

## get rid of poor temporal resolution
pbdb.dat <- pbdb.dat[pbdb.dat$collections.10_my_bin != "",]

##	get rid of bad taxonomy
pbdb.dat <- pbdb.dat[pbdb.dat$occurrences.order_name != "",]
pbdb.dat <- pbdb.dat[-which(pbdb.dat$occurrences.order_name == 'Ammonitida'),]

##	drop missing levels
pbdb.dat$collections.10_my_bin <- as.factor(as.character(pbdb.dat$collections.10_my_bin))
pbdb.dat$occurrences.order_name <- as.factor(as.character(pbdb.dat$occurrences.order_name))
pbdb.dat$occurrences.genus_name <- as.factor(as.character(pbdb.dat$occurrences.genus_name))
pbdb.dat$collections.reference_no <- as.factor(as.character(pbdb.dat$collections.reference_no))
pbdb.dat$collection_no <- as.factor(as.character(pbdb.dat$collection_no))


##	subsampled diversity (for comparison's sake)
pbdb.samp <- read.csv("subsampled_curve_data.csv")

##	raw diversity curve (for 3 timer stat, etc)
pbdb.curv <- read.csv("raw_curve_data.csv")

##	get bin times
pbdb.time <- pbdb.samp$Midpoint.Ma
names(pbdb.time) <- pbdb.samp$Bin.name
pbdb.time <- pbdb.time[levels(pbdb.samp$Bin.name)]

##  data.frame of publication, diversity and 3T stat
ord.tbin.bias <- aggregate(list(div=pbdb.dat$occurrences.genus_name),
						   list(ord=pbdb.dat$occurrences.order_name,
						   		tbin=pbdb.dat$collections.10_my_bin),
						   function(x) length(unique(x)))

ord.tbin.bias$T3.stat <- pbdb.curv$Three.timer.sampling.stat[match(ord.tbin.bias$tbin,pbdb.curv$Bin.name)]
ord.tbin.bias$T3.div <- ord.tbin.bias$div/ord.tbin.bias$T3.stat
head(ord.tbin.bias)

##	record pubs per tbin
tbin.pub <- tapply(pbdb.dat$collections.reference_no,pbdb.dat$collections.10_my_bin,function(x) length(unique(x)))
ord.tbin.bias$tbin.pub <- tbin.pub[ord.tbin.bias$tbin]

######### PLOTTING!!!! ############
setwd('../../../code/')
## makes `figSupp_divByPubOrd.pdf'

##	calculate corrected diversity
if(makePlot) {
	source('sstat_plotting_par.R')
	with(plot.pars, {
		quartz(width=width,height=height)
		par(mar=mar,mgp=mgp)
	})
}
pbdb.ord.div <- with(ord.tbin.bias,
	pbdb.3t.pub(div, T3.stat, tbin.pub, ord, tbin, pbdb.time, min.pub=10, plotit=makePlot)
)
rownames(pbdb.ord.div)

##	function `pbdb.3t.pub' modified to save regression to global env
pbdb.pub.lm

##	using correction on genus occurances
pbdb.gen.occ <- with(pbdb.dat,samp2site.spp(collections.10_my_bin,occurrences.genus_name,rep(1,nrow(pbdb.dat))))
pbdb.gen.occ <- pbdb.gen.occ[rownames(pbdb.ord.div),]
pbdb.gen.occ <- pbdb.gen.occ[,colSums(pbdb.gen.occ) > 0]
pbdb.gen.occ <- 1*(pbdb.gen.occ > 0)

# # fill-in Lazzarous taxa...might not be best idea
# pbdb.gen.occ <- apply(pbdb.gen.occ,2,function(x) {
	# occInd <- which(x > 0)
	# x[min(occInd):max(occInd)] <- 1
	# return(x)
# })

##  data.frame for predicting from pbdb.pub.lm
pub.data <- with(ord.tbin.bias,data.frame(log(tbin.pub[match(rownames(pbdb.gen.occ),tbin)])))
rownames(pub.data) <- NULL
colnames(pub.data) <- names(pbdb.pub.lm$coeff)[2]

##  corrected genus-level data
pbdb.gen.occ3TP <- pbdb.gen.occ / pbdb.curv$Three.timer.sampling.stat[match(rownames(pbdb.gen.occ),pbdb.curv$Bin.name)]
pbdb.gen.occ3TP <- pbdb.gen.occ3TP * exp(-predict(pbdb.pub.lm,newdata=pub.data))

# gen2ord <- as.character(with(pbdb.dat,
	# occurrences.order_name[match(colnames(pbdb.gen.occ3TP),occurrences.genus_name)]))


##  data.frame of different diversity measures
pbdb.div <- aggregate(list(raw=pbdb.dat$occurrences.genus_name),
					   list(tbin=pbdb.dat$collections.10_my_bin),
					   function(x) length(unique(x)))
pbdb.div$time <- pbdb.time[as.character(pbdb.div$tbin)]

pbdb.div$sqs <- pbdb.samp[match(pbdb.div$tbin,pbdb.samp$Bin.name), "Mean.sampled.diversity"]

pbdb.div$pub3t <- rowSums(pbdb.ord.div)[match(pbdb.div$tbin,rownames(pbdb.ord.div))]

pbdb.div <- pbdb.div[order(pbdb.div$time, decreasing=TRUE),]

########### PLOTTING!!!!! ##########
# makes `figSupp_sqsRaw3tpub.pdf'
if(makePlot) {
	source('sstat_plotting_par.R')
	with(plot.pars, {
		# quartz(width=width*2.1,height=height)
		par(mar=mar, mgp=mgp, mfrow=c(1,2))

		with(pbdb.div[2:48,], {
			paleoPlot(time[-1], scale(diff(pub3t)),type='l',col='red',
					  y.lim=c(-3,3),ylab="Scaled diversity fluctuations")
			lines(time[-1], scale(diff(sqs)))
			lines(time[-1], scale(diff(raw)), lty=2)
			text(par("usr")[2],par("usr")[4],labels="A",adj=adj,cex=cex.txt)
			
			logPlot(my.ecdf(abs(scale(pub3t)),complement=TRUE),type="l",log="xy",col="red",
					xlab="|Scaled fluctuations|", ylab="Cumulative density")
			lines(my.ecdf(abs(scale(sqs)),complement=TRUE))
			lines(my.ecdf(abs(scale(raw)),complement=TRUE),lty=2)
			text(10^par("usr")[2],10^par("usr")[4],labels="B",adj=adj,cex=cex.txt)
		})
	})
}


##	class-level diversity
# match orders to classes
ord2cls <- with(pbdb.dat, occurrences.class_name[match(colnames(pbdb.ord.div),occurrences.order_name)])
ord2cls <- as.character(ord2cls)
names(ord2cls) <- colnames(pbdb.ord.div)  # 20 orders not classified to class

# aggregate diversity for classes
pbdb.cls.div <- sapply(split(data.frame(t(pbdb.ord.div)),ord2cls), colSums)


##	corrected flux
pbdb.ord.flux <- apply(pbdb.ord.div,2,function(x)
							{
								raw.flux <- diff(c(0,x))
								return(raw.flux[raw.flux != 0])
							})
# (same just log-transform div first)
pbdb.ord.flux.log <- apply(pbdb.ord.div,2,function(x)
							{
								raw.flux <- diff(log(x))
								bad <- is.na(raw.flux) | !is.finite(raw.flux)
								return(raw.flux[!bad])
							})


pbdb.cls.flux <- apply(pbdb.cls.div,2,function(x)
							{
								raw.flux <- diff(c(0,x))
								return(raw.flux[raw.flux != 0])
							})

##	sstat analysis
pbdb.sstat.ord.cor <- sstat.comp(pbdb.ord.flux,minN=10,plotit=makePlot)
pbdb.sstat.ord.log <- sstat.comp(pbdb.ord.flux.log,minN=10,plotit=makePlot)

if(makePlot) {
plot(my.ecdf(scale(log(pbdb.sstat.ord.log$beta))))
curve(pnorm(x),add=TRUE)

library(MASS)
plot(my.ecdf(pbdb.sstat.ord.log$beta),log='x')
with(pbdb.sstat.ord.log, curve(pgamma(x,gam.par[1],gam.par[2]),add=TRUE,col='red'))

lnorm.par <- with(pbdb.sstat.ord.log, fitdistr(beta,'lognormal'))
with(lnorm.par, curve(plnorm(x,estimate[1],estimate[2]),add=TRUE,col='blue'))

chisq.par <- with(pbdb.sstat.ord.log, fitdistr(beta,'chi-squared',start=list(df=1,ncp=1))) 
with(chisq.par, curve(pchisq(x,estimate[1],estimate[2]),add=TRUE,col='goldenrod3'))

Px.lnorm.prm <- function(b,x,logm,logsd) {
	dnorm(x,sd=b^-0.5)*dlnorm(b,logm,logsd)
}

x <- c(seq(-20,-4,length=100),seq(-4,4,length=500),seq(4,20,length=100))
Px.lnorm <- sapply(seq(-20,20,length=500), function(x) {
	with(lnorm.par,integrate(Vectorize(Px.lnorm.prm),lower=0,upper=Inf,x=x,logm=estimate[1],logsd=estimate[2]))$value
})

}

# plot(my.ecdf(pbdb.sstat.ord.cor$beta,complement=TRUE),log="x")
# with(pbdb.sstat.ord.cor, curve(pgamma(x,gam.par[1],gam.par[2],lower.tail=FALSE),add=TRUE))

# plot(pbdb.sstat.ord.cor)

# plot(do.call(rbind,lapply(pbdb.sstat.ord.cor$Px.sub,function(x) my.ecdf(scale(x)[,1],TRUE))))
# curve(pnorm(x,lower.tail=FALSE),col="red",add=TRUE)

pbdb.sstat.cls.cor <- sstat.comp(pbdb.cls.flux,minN=10,plotit=makePlot)

setwd(oldwd)




