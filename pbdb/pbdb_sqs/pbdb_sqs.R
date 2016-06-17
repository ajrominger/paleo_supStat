library(vegan)

pbdb.sqs <- function(dat,q,pars,B=50,raw.id=FALSE) {
	if(missing(pars)) pars <- calc.sqs.par(dat)
	
	if(raw.id) {
		out <- vector("list",B)
		
		for(i in 1:B) {
			out[[i]] <- calc.sqs(dat,q,pars,raw.id=raw.id)
		}
		
		return(list(ID=lapply(out,function(x) x$colID),div.data=out))
	} else {
		out <- matrix(NA,nrow=3,ncol=B)
		
		for(i in 1:B) {
			out[i,] <- calc.sqs(dat,q,pars,raw.id=raw.id)
		}
		
		##	summary on out (rows of out are div, freq, new.q)
		
		return(apply(out,2,function(x) c(mean(x), quantile(x,c(0.025,0.975)))))
	}
}

calc.sqs.par <- function(dat) {
	## taxa.occ is vector of all occurances by taxon combined across pubs and colls
	taxa.occ <- as.numeric(table(dat$occurrences.genus_name))
	
	## taxa.pub is a vector of count of pubs for each taxon
	taxa.pub <- tapply(dat$collections.reference_no,dat$occurrences.genus_name,function(x) length(unique(x)))
	
	## colXtaxa is a matrix of presence/absence with collections=rows and taxa=columns
	colXtaxa <- do.call(rbind,tapply(as.factor(dat$occurrences.genus_name),dat$collection_no,table,simplify=FALSE))
	colXtaxa[colXtaxa > 0] <- 1		# just to be sure that matrix is presence/absence
	
	## taxa.col is count of collections for each taxon
	taxa.col <- apply(colXtaxa,2,sum)
	## most.div is vector of presence/absence in most diverse collection
	most.div <- colXtaxa[which.max(apply(colXtaxa,1,sum)),]
	
	O <- sum(taxa.occ)		# O in alroy 2010 science
	n.max <- max(taxa.occ)	# n_1 in alroy 2010 science
	p1 <- sum(taxa.occ[taxa.pub == 1])	# p1 in alroy 2010 science
	
	# t_max in alroy 2010 science
	n.div <- sum(taxa.occ[most.div == 1 & taxa.col == 1])
	
	# u' in alroy 2010 science
	u.prm <- (O - n.max - p1 + n.div) / (O - n.max)
	if(!is.finite(u.prm)) {
		u.prm <- 0	# happens if denominator = 0 & numerator = 0 (e.g. 1 occurence from 1 collection and 1 pub)
	} else if(u.prm < 0) {
		u.prm <- 0	# happens if (n.max + p1) > (O + n.div) 
	}
	
	return(list(O=O, n.max=n.max, taxa.col=taxa.col, colXtaxa=colXtaxa, p1=p1,
				taxa.occ=taxa.occ, taxa.pub=taxa.pub, n.div=n.div, u.prm=u.prm))
}

calc.sqs <- function(dat,q,pars,raw.id=FALSE) {
#	time.sofar <- proc.time()[3]
	
#	attach(pars,pos=-1)
	u.prm <- pars$u.prm
	O <- pars$O
	n.max <- pars$n.max
	colXtaxa <- pars$colXtaxa
	taxa.occ <- pars$taxa.occ
	
	
	if(q > u.prm) {
		cat("q set higher than coverage as estimated by Good's u with Alroy's puclication correction")
		return(c(NA,NA,NA))
	}
	
	new.q <- q / u.prm	# smaller the coverage, harder the sampling (i.e. bigger q)
	
	freq <- taxa.occ / (O - n.max)	# shares of each taxon (ignoring `abund' of most common)
	col.taxa.share <- t(t(colXtaxa)*freq)	# matrix of collections by taxa with cells = each taxon's share
	
#	print(proc.time()[3] - time.sofar)
#	time.sofar <- proc.time()[3]
	
	##	estimate of how shares accumulate with subsampled collections
	share.accum <- specaccum(colXtaxa,method="random",permutations=1)$richness*mean(freq)
	expt.quorum <- max(which(share.accum - new.q <= 0))    # expected number of collections hitting quorum
	
#	print(proc.time()[3] - time.sofar)
#	time.sofar <- proc.time()[3]
	
	##	begin sampling
	nCol <- nrow(colXtaxa)	# number of collections
	col.samp.id <- sample(nCol,expt.quorum)	# initial try at which collections to sample
	
	## calulate sampled frequency/coverage
	samp.freq <- apply(col.taxa.share[col.samp.id,],2,sum)/apply(colXtaxa[col.samp.id,],2,sum)
	samp.freq[!is.finite(samp.freq)] <- 0
	samp.freq <- sum(samp.freq)
	
	##	if too many, cut number of collections by proportion new.q/samp.freq
	if(samp.freq > new.q) {
		rm.these <- sample(length(col.samp.id),size=ceiling(expt.quorum*new.q/samp.freq))
		col.samp.id <- col.samp.id[-rm.these]
		
		samp.freq <- apply(col.taxa.share[col.samp.id,],2,sum)/apply(colXtaxa[col.samp.id,],2,sum)
		samp.freq[!is.finite(samp.freq)] <- 0
		samp.freq <- sum(samp.freq)
	}
	
	print(samp.freq)
	
#	print(proc.time()[3] - time.sofar)
#	time.sofar <- proc.time()[3]
	
	if(samp.freq < new.q) {			# add collections
		print("too few")
		i <- 1
		while(samp.freq < new.q) {
			print(i)
			## add one new collection
			col.samp.id <- c(col.samp.id,
							 sample((1:nCol)[!(1:nCol %in% col.samp.id)],1))
			samp.freq <- apply(col.taxa.share[col.samp.id,],2,sum)/apply(colXtaxa[col.samp.id,],2,sum)
			
			samp.freq[!is.finite(samp.freq)] <- 0
			samp.freq <- sum(samp.freq)
			
			i <- i + 1
		}
	} else if(samp.freq > new.q) {	# take away collections
		print("too many")
		i <- 1
		while(samp.freq > new.q) {
			print(i)
			## remove one collection
			col.samp.id <- col.samp.id[-sample(length(col.samp.id),1)]
			samp.freq <- apply(col.taxa.share[col.samp.id,],2,sum)/apply(colXtaxa[col.samp.id,],2,sum)
			
			samp.freq[!is.finite(samp.freq)] <- 0
			samp.freq <- sum(samp.freq)
			
			i <- i + 1
		}
		
		##	now sampled coverage is lower, so add back to be consistent
		i <- 1
		while(samp.freq < new.q) {
			print(i)
			## add one new collection
			col.samp.id <- c(col.samp.id,
							 sample((1:nCol)[!(1:nCol %in% col.samp.id)],1))
			samp.freq <- apply(col.taxa.share[col.samp.id,],2,sum)/apply(colXtaxa[col.samp.id,],2,sum)
			
			samp.freq[!is.finite(samp.freq)] <- 0
			samp.freq <- sum(samp.freq)
			
			i <- i + 1
		}
	}
	
	##	total diversity sampled
	samp.div <- sum(apply(colXtaxa[col.samp.id,],2,sum) > 0)
	
#	detach(pars,pos=-1)
	if(raw.id) {
		return(list(div=samp.div,coverage=samp.freq,new.q=new.q,colID=as.numeric(rownames(colXtaxa)[col.samp.id])))
	} else {
		return(c(div=samp.div,coverage=samp.freq,new.q=new.q))
	}
}

