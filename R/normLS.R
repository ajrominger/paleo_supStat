normLS <- function(data,comp=FALSE) {
	par.init <- c(mean(data),sd(data))
	optim(par.init,norm.ss,data=data,comp=comp)
}

norm.ss <- function(pars,data,comp) {
	mean <- pars[1]
	sd <- pars[2]
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
	difz <- pnorm(xval,mean=mean,sd=sd,lower.tail=lower) - yval
	sum(difz^2)
}
