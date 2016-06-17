#########  calculate coverage of pbdb data and play with sqs  #########

##	load data
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb_read-in.R")

##	load functions
source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb_sqs.R")

##	some things easier as char (otherwise will need to drop levels in apply below)
pbdb.dat$occurrences.genus_name <- as.character(pbdb.dat$occurrences.genus_name)
pbdb.dat$collection_no <- as.character(pbdb.dat$collection_no)
pbdb.dat$collections.reference_no <- as.character(pbdb.dat$collections.reference_no)

##	extract indices of each order-10mybin combo
pbdb.ord.tbin <- aggregate(1:nrow(pbdb.dat),list(ord=as.character(pbdb.dat$occurrences.order_name),
												 tbin=as.character(pbdb.dat$collections.10mybin)),c)
apply(pbdb.ord.tbin,1,function(X) length(unique(pbdb.dat[X$x,"collection_no"])))


##	compute coverage
pbdb.ord.sqs.par <- apply(pbdb.ord.tbin,1,function(X)
							{
								dat <- pbdb.dat[X$x,]
								calc.sqs.par(dat)
							})

pbdb.ordT.u <- sapply(pbdb.ord.sqs.par,function(x) x$u.prm)
pbdb.ordT.n <- sapply(pbdb.ord.tbin[,3],length)
plot(pbdb.ordT.n,pbdb.ordT.u,log="x")
abline(v=25)

hist(pbdb.ordT.u,breaks=c(0,0.01,seq(0.1,1,by=0.1)),xlab="Good's u (i.e. coverage)",main="",prob=FALSE)

barplot(table(pbdb.ordT.u >= 0.5),names.arg=c("< 0.5",">= 0.5"),space=0.2,xlab="Good's u")


##	matrix of time X order with cells = u, and cells = n (occurrences)
tbinXord.u <- matrix(0,nrow=length(unique(pbdb.ord.tbin[,2])),
					  ncol=length(unique(pbdb.ord.tbin[,1])))

rownames(tbinXord.u) <- unique(pbdb.ord.tbin[,2])
colnames(tbinXord.u) <- unique(pbdb.ord.tbin[,1])

tbinXord.n <- tbinXord.u

##	probably there is a better way to do this...
for(i in 1:nrow(pbdb.ord.tbin)) {
	tbinXord.u[pbdb.ord.tbin[i,2],pbdb.ord.tbin[i,1]] <- pbdb.ordT.u[i]
	tbinXord.n[pbdb.ord.tbin[i,2],pbdb.ord.tbin[i,1]] <- pbdb.ordT.n[i]
}


##	matrix of T/F whether an order (col) in a time period (row) has u >= u.crit & n >= n.crit
u.crit <- 0.5
n.crit <- 20
tbinXord.use <- tbinXord.u >= u.crit & tbinXord.n >= n.crit

##	number of orders w/ u >= u.crit & n >= n.crit & num time bins >= ntime.crit
ntime.crit <- 10
#hist(apply(tbinXord.use,2,sum))
ord2use <- apply(tbinXord.use,2,sum) >= ntime.crit
sum(ord2use)

logPlot(tbinXord.n,tbinXord.u,log="x",xlab="Number of occurrences",ylab="Good's u")
points(tbinXord.n[,ord2use][tbinXord.use[,ord2use]],tbinXord.u[,ord2use][tbinXord.use[,ord2use]],col="red")



pbdb.ord.nU <- apply(tbinXord.u,2,function(x) sum(x >= 0.5))
hist(pbdb.ord.nU,xlab="Number of time periods per order\nwith u >= 0.5",main="")

barplot(table(pbdb.ord.nU >= 10),names.arg=c("< 10",">= 10"),space=0.2,xlab="Number of time periods per order\nwith u >= 0.5")

pbdb.ord.nN <- apply(tbinXord.n,2,function(x) sum(x >= 20))

##	finally, make object like pbdb.ord.tbin, but only having those orders with >= 10 points at u >= 0.5
these.tbin10 <- which(pbdb.ordT.u >= 0.5 & pbdb.ordT.n >= 20)
pbdb.ord.tbin[these.tbin10,1] 

pbdb.ord.tbin10 <- pbdb.ord.tbin[pbdb.ordT.u >= 0.5 & pbdb.ordT.n >= 20,]   # get only rows associated with u >= 0.5 & n >= 20
# now take only orders with >= 10 points
pbdb.ord.tbin10 <- pbdb.ord.tbin10[pbdb.ord.tbin10[,1] %in% colnames(tbinXord.u)[pbdb.ord.nU >= 10],]

nrow(pbdb.ord.tbin10)
sapply(pbdb.ord.tbin10[,3],length)

##	do same for sqs pars
pbdb.ord.sqs.par10 <- pbdb.ord.sqs.par[pbdb.ordT.u >= 0.5 & pbdb.ordT.n >= 20]