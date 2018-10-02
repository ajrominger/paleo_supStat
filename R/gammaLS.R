gammaLS <- function(data,comp=FALSE) {
	par.init <- c(mean(data)^2,mean(data))/var(data)
	optim(par.init,gamma.ss,data=data,comp=comp)
}

gamma.ss <- function(pars,data,comp) {
	shape <- pars[1]
	rate <- pars[2]
	tabz <- table(data)
	xval <- as.numeric(names(tabz))
	yval <- cumsum(as.numeric(tabz))/sum(tabz)
	if(comp) {
		yval <- 1 - yval
		yval <- c(1,yval[-length(yval)])
		lower <- FALSE
	} else {
		lower <- TRUE
	}
	difz <- pgamma(xval,shape=shape,rate=rate,lower.tail=lower) - yval
	sum(difz^2)
}
