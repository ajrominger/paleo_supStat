##############  make P(x) fig for PBDB raw to go in supplament  ##############

##	load data
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb/pbdb_read-in.R")

##	load functions
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_comp.R")
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/code/sstat_methods.R")


##	directory to save figs
save.dir <- "/Users/andrewrominger/Desktop/Research/paleo_supStat/ms/current_version"

pbdb.ord.sstat <- sstat.comp(pbdb.ord.flux,minN=15,plotit=FALSE)

##	plotting
quartz(width=4,height=4)
par(mar=c(3.5,3.5,1,0.5)+0.1,mgp=c(2.5,1,0))
plot(pbdb.ord.sstat,xlab="|Fluctuations|",ylab="Cumulative density")

dev.print(pdf,file=paste(save.dir,"figSupp_pbdbRaw_Px.pdf",sep="/"),width=4,height=4)