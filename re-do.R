#######  corrects diversity by 3-timer then number of pubs  #######
makePlot <- TRUE


oldwd <- setwd('~/Dropbox/Research/paleo_supStat')

##  convenience function to produce a matrix of time by ord with cells
##  of corrected diversity
source('code/pbdb_3t_pub.R')

##	load other needed funcitons
paleoPlot <- function(...) {
    plot(..., xaxt = 'n')
    socorro::paleoAxis(1)
}
# source("~/R_functions/paleoPlot.R")

samp2site.spp <- socorro::tidy2mat
# source("~/R_functions/samp2site_spp.R")

logPlot <- function(..., log) {
    plot(..., log = log, 
         xaxt = ifelse(grepl('x', log), 'n', 's'), 
         yaxt = ifelse(grepl('y', log), 'n', 's'))
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

##	record pubs per tbin
tbin.pub <- tapply(pbdb.dat$collections.reference_no,pbdb.dat$collections.10_my_bin,function(x) length(unique(x)))
ord.tbin.bias$tbin.pub <- tbin.pub[ord.tbin.bias$tbin]

######### PLOTTING!!!! ############
setwd('../../../code/')
## makes `figSupp_divByPubOrd.pdf'

##	calculate corrected diversity

pbdb.ord.div <- with(ord.tbin.bias,
                     pbdb.3t.pub(div, T3.stat, tbin.pub, ord, tbin, pbdb.time, min.pub=10, plotit=makePlot))

##	function `pbdb.3t.pub' modified to save regression to global env
pbdb.pub.lm

##	using correction on genus occurances
pbdb.gen.occ <- with(pbdb.dat,samp2site.spp(collections.10_my_bin,occurrences.genus_name,rep(1,nrow(pbdb.dat))))
pbdb.gen.occ <- pbdb.gen.occ[rownames(pbdb.ord.div),]
pbdb.gen.occ <- pbdb.gen.occ[,colSums(pbdb.gen.occ) > 0]
pbdb.gen.occ <- 1*(pbdb.gen.occ > 0)

##  data.frame for predicting from pbdb.pub.lm
pub.data <- with(ord.tbin.bias,data.frame(log(tbin.pub[match(rownames(pbdb.gen.occ),tbin)])))
rownames(pub.data) <- NULL
colnames(pub.data) <- names(pbdb.pub.lm$coeff)[2]

##  corrected genus-level data
pbdb.gen.occ3TP <- pbdb.gen.occ / pbdb.curv$Three.timer.sampling.stat[match(rownames(pbdb.gen.occ),pbdb.curv$Bin.name)]
pbdb.gen.occ3TP <- pbdb.gen.occ3TP * exp(-predict(pbdb.pub.lm,newdata=pub.data))


##  data.frame of different diversity measures
pbdb.div <- aggregate(list(raw=pbdb.dat$occurrences.genus_name),
                      list(tbin=pbdb.dat$collections.10_my_bin),
                      function(x) length(unique(x)))
pbdb.div$time <- pbdb.time[as.character(pbdb.div$tbin)]

pbdb.div$sqs <- pbdb.samp[match(pbdb.div$tbin,pbdb.samp$Bin.name), "Mean.sampled.diversity"]

pbdb.div$pub3t <- rowSums(pbdb.ord.div)[match(pbdb.div$tbin,rownames(pbdb.ord.div))]

pbdb.div <- pbdb.div[order(pbdb.div$time, decreasing=TRUE),]

##	corrected flux
pbdb.ord.flux <- apply(pbdb.ord.div,2,function(x)
{
    raw.flux <- diff(c(0,x))
    return(raw.flux[raw.flux != 0])
})

##	sstat analysis
pbdb.sstat.ord.cor <- sstat.comp(pbdb.ord.flux,minN=10,plotit=makePlot)
