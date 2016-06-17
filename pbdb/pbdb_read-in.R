pbdb.dat <- data.frame()

for(i in 1:4) {
	pbdb.dat <- rbind(pbdb.dat,read.csv(paste("~/Desktop/Fulbright/paleo_superStat/pbdb/pbdb_data/metazoa",i,"-occs.csv",sep="")))
}

##	missing orders to NA
pbdb.dat$occurrences.order_name[pbdb.dat$occurrences.order_name == ""] <- NA
pbdb.dat$occurrences.order_name <- as.factor(as.character(pbdb.dat$occurrences.order_name))  # drop NA from levels

##	poor temporal resolution deleted
pbdb.dat <- pbdb.dat[pbdb.dat$collections.10mybin != "",]
pbdb.dat$collections.10mybin <- as.factor(as.character(pbdb.dat$collections.10mybin))  # drop "" from levels

##	drop missing levels in orders and genera
pbdb.dat$occurrences.order_name <- as.factor(as.character(pbdb.dat$occurrences.order_name))
pbdb.dat$occurrences.genus_name <- as.factor(as.character(pbdb.dat$occurrences.genus_name))

##	calculate diversity per order, per time
pbdb.ord.div <-  aggregate(list(gen=pbdb.dat$occurrences.genus_name),
						   list(ord=pbdb.dat$occurrences.order_name,tbin=pbdb.dat$collections.10mybin),
						   function(x) length(unique(x)))

##	get measure of times for each bin (really should look this up in PBDB)
pbdb.time <- tapply(pbdb.dat$collections.ma_mid,pbdb.dat$collections.10mybin,mean)

##	tag time onto pbdb.ord.div object
pbdb.ord.div$time <- pbdb.time[as.character(pbdb.ord.div$tbin)]
pbdb.ord.div <- pbdb.ord.div[order(pbdb.ord.div$time,decreasing=TRUE),]



pbdb.ord.divflux <- split(pbdb.ord.div,pbdb.ord.div$ord,drop=FALSE)
pbdb.ord.divflux <- lapply(pbdb.ord.divflux,function(x)
						{
							x$flux <- c(x$gen[1],diff(x$gen))
							x$flux.cent <- x$flux - mean(x$flux)
							colnames(x) <- c("ord","tbin","div","time","flux","flux.cent")
							return(x)
						})

pbdb.ord.divflux2 <- lapply(pbdb.ord.divflux,function(x)
							{
								newx <- x[,c("time","div","flux")]
								out <- as.data.frame(cbind(time=sort(pbdb.time,decreasing=TRUE),div=0,flux=0))
								out[match(x$tbin,rownames(out)),2:3] <- newx[,2:3]
								
								return(out)
							})


pbdb.ord.divMa <- matrix(0,length(pbdb.time),length(pbdb.ord.divflux))
rownames(pbdb.ord.divMa) <- names(pbdb.time)
colnames(pbdb.ord.divMa) <- names(pbdb.ord.divflux)

for(i in 1:length(pbdb.ord.divflux)) {
	pbdb.ord.divMa[pbdb.ord.divflux[[i]]$tbin,i] <- pbdb.ord.divflux[[i]]$div
}

pbdb.ord.divMa <- pbdb.ord.divMa[,apply(pbdb.ord.divMa,2,function(x) sum(x > 0)) >= 10]
write.table(pbdb.ord.divMa,file="~/Desktop/pbdb_ord_div.txt",sep="\t",row.names=FALSE)

##	note: this is taking the re-centered data
pbdb.ord.flux <- lapply(pbdb.ord.divflux,function(x) x$flux.cent[x$flux != 0])


