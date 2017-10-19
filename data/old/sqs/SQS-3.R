
# sqs version 3.0 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# can be used to perform classical rarefaction instead of SQS
# written 29 July 2010; version 2.0 completed 14 February 2011; version 3.0 written
#  3 June 2011
# changes in version 3.0: an even better subsampling algorithm involving a new
#  adjustment to u combined with a new throwback criterion
# changes in version 2.0: improved subsampling algorithm; including the dominant 
#  taxon is now the default; improved reporting of errors and basic statistics
# warning: do not use this program with taxonomic occurrence data drawn from
#  multiple published references because it is not designed to count
#  single-reference taxa or adjust for long taxonomic lists
# warning: version 1.0 yields estimates that are downwards-biased when q < 0.6
#  and abundance distributions are highly uneven

sqs<-function(ab,q,trials=100,ex.dom=TRUE) {

	params <- array(data=NA,dim=0,dimnames=c("raw richness"))

	# compute basic statistics
	specimens <- sum(ab)		# total number of speciments
	singletons <- sum(ab == 1)	# number of singletons
	doubletons <- sum(ab == 2)	# number of doubletons
	highest <- max(ab)			# max frequency
	mostfrequent <- which(ab == highest)	# most frequent taxon
	
	
	# calculate u, excluding or including dominate taxon
	if (ex.dom) {
		u <- 1 - singletons / (specimens - highest)  # should always be >= 0
	} else {
		u <- 1 - singletons / specimens
	}

	if (u == 0)	{
		stop("Coverage is zero because all taxa are singletons")
	}

	# compute raw taxon frequencies (temporarily)
	freq <- ab / specimens

	params["raw richness"] <- length(ab)
	params["Good's u"] <- u
	params["subsampled richness"] <- NA
	params["subsampled u"] <- NA


	params["dominance"] <- highest / specimens
	params["specimens"] <- specimens
	params["singletons"] <- singletons
	params["doubletons"] <- doubletons
	params["specimens drawn"] <- 0

	if (!ex.dom) {
		highest <- 0
		mostfrequent <- 0
	}

	# return if the quorum target is higher than overall coverage
	if (q > u) {
		cat("quorum set higher than coverate")
		return(params)
	}

	# compute raw taxon frequencies...should there be an option to not have '- highest'?
	freq <- ab / (specimens - highest)

	# create an array in which each cell corresponds to one specimen
	ids <- rep(1:length(ab),ab)

	# subsampling trial loop
	# s will be the subsampled taxon count
	s <- array(rep(0,trials))
	subsingle <- array(rep(0,trials))
	subdouble <- array(rep(0,trials))
	subchao <- array(rep(0,trials))
	mostfrequentdrawn <- 0

	for(trial in 1:trials)	{
		pool <- ids
		left <- length(pool)
		seen <- array(data=rep(0,length(ab)))
		subfreq <- array(rep(0,length(ab)))
		if(method != "rarefaction" && method != "CR") {
			udrawn <- 0
			# equation new to version 3.0
			# the exponent corrects for downwards bias
			while (udrawn**u * u < q)	{
			# draw a specimen
				x <- floor(runif(1,min=1,max=left+1))
			# add to frequency and taxon sums if species has
			#  not been drawn previously
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if(seen[pool[x]] == 0)	{
					if(pool[x] != mostfrequent) {
						udrawn <- udrawn + freq[pool[x]]
					}
					seen[pool[x]] <- 1
			# randomly throw back some draws that put the sum over q
			#  (an even better algorithm added in version 3.0)
					if (runif(1) <= freq[pool[x]] || udrawn**u * u < q) {
						s[trial] <- s[trial] + 1
					} else {
						subfreq[pool[x]] <- subfreq[pool[x]] - 1
					}
				}
			# decrease pool of specimens not yet drawn
				pool[x] <- pool[left]
				left <- left - 1
			}
		} else {
			i <- 0
			draws <- 0
			while (i < q)	{
				draws <- draws + 1
				x <- floor(runif(1,min=1,max=length(ids)-draws+2))
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (pool[x] != mostfrequent)	{
					i <- i + 1
				}
				if (seen[pool[x]] == 0)	{
					seen[pool[x]] <- 1
					s[trial] <- s[trial] + 1
				}
				pool[x] <- pool[length(ids)-draws+1]
			}
		}
		for(i in 1:length(ab)) {
			if(subfreq[i] == 1 && i != mostfrequent) {
				subsingle[trial] <- subsingle[trial] + 1
			} else if(subfreq[i] == 2 && i != mostfrequent) {
				subdouble[trial] <- subdouble[trial] + 1
			}
		}
		if (subsingle[trial] > 0 && subdouble[trial] > 0) {
			subchao[trial] <- s[trial] + subsingle[trial]**2/(2*subdouble[trial])
		} else {
			subchao[trial] <- s[trial]
		}
		params["specimens drawn"] <- params["specimens drawn"] + sum(subfreq)
		if(mostfrequent != 0) {
			mostfrequentdrawn <- mostfrequentdrawn + subfreq[mostfrequent]
		}
	}

	params["specimens drawn"] <- params["specimens drawn"] / trials
	# compute vector of non-zero counts
	options(warn=-1)
	s2 <- sort(sqrt(s-1))^2+1
	options(warn=0)
	# compute geometric mean
	params["subsampled richness"] <- exp(mean(log(s2))) * length(s2)/length(s)
	# use of arithmetic means to compute Good's u is adequate
	mostfrequentdrawn <- mostfrequentdrawn / trials
	params["subsampled u"] <- 1 - mean(subsingle) / (params["specimens drawn"] - mostfrequentdrawn)
	return(params)

}

