#############  permatuation test of order-level uniqueness (PBDB)  #############

##	convenience funcitons
source("~/R_functions/ks_stat_pfun.R")
source("~/R_functions/den_fill.R")

##	create matrix of taxa (gen) by time
make.occ.ma <- function(x) {	# conenience function to creat matrix
	x$occurrences.genus_name <- as.factor(as.character(x$occurrences.genus_name))
	y <- tapply(rep(1,nrow(x)), list(x$occurrences.genus_name, x$collections.10mybin), function(X) as.numeric(sum(X) > 0))
	y[is.na(y)] <- 0
	
	return(as.data.frame(y))
}

##	extract fluctuations (specific for PBDB)
make.pbdb.flux <- function(occ.ma,ord) {
	ord.split <- split(occ.ma,ord)
	
	lapply(ord.split,function(x)
				{
					flux <- diff(c(0,apply(x,2,sum)))
					flux <- flux[flux != 0]
					return(flux)
				})
}


##	depends on .../pbdb_readIn.R

##	get rid of missing orders
pbdb.dat2 <- pbdb.dat[!is.na(pbdb.dat$occurrences.order_name),]
pbdb.dat2$occurrences.order_name <- as.factor(as.character(pbdb.dat2$occurrences.order_name))
pbdb.dat2$occurrences.genus_name <- as.factor(as.character(pbdb.dat2$occurrences.genus_name))


##	matrix of occurences
pbdb.occ.ma <- make.occ.ma(pbdb.dat2)
pbdb.occ.ma <- pbdb.occ.ma[,order(pbdb.time,decreasing=TRUE)]

# should be TRUE
all(rownames(pbdb.occ.ma) %in% as.character(pbdb.dat2$occurrences.genus_name))

pbdb.occ.ord <- pbdb.dat2$occurrences.order_name[match(rownames(pbdb.occ.ma),pbdb.dat2$occurrences.genus_name)]

# should be FALSE
any(is.na(pbdb.occ.ord))


##	calculate real sstat object
pbdb.ord.flux <- make.pbdb.flux(pbdb.occ.ma,pbdb.occ.ord)
pbdb.ord.sstat <- sstat.comp(pbdb.ord.flux,minN=15)
pbdb.raw.PPx <- function(x) 0.5 + 0.5*pbdb.ord.sstat$PPx(x,FALSE)
pbdb.ord.D <- ks.stat.pfun(unlist(pbdb.ord.sstat$raw.pk),"pbdb.raw.PPx")

##	begin permutation

nrun <- 100
pbdb.ord.permD <- numeric(nrun)
pbdb.ord.perm.sstat <- vector("list",nrun)

for(i in 1:nrun) {
	##	permute orders
	this.occ.ord <- sample(pbdb.occ.ord)
	
	##	re-calculate flux by order
	this.ord.flux <- make.pbdb.flux(pbdb.occ.ma,this.occ.ord)
	
	##	calc super-stat and store it
	this.sstat <- sstat.comp(this.ord.flux,minN=15,plotit=FALSE)
	pbdb.ord.perm.sstat[[i]] <- this.sstat
	
	##	calc D-statistic and store it
	this.PPx <- function(x) 0.5 + 0.5*this.sstat$PPx(x,FALSE)
	
	pbdb.ord.permD[i] <- ks.stat.pfun(unlist(this.sstat$raw.pk),"this.PPx")
	
	print(i)
}


den.fill(pbdb.ord.permD,xlim=range(pbdb.ord.permD,pbdb.ord.D))
abline(v=pbdb.ord.D,col="red")
