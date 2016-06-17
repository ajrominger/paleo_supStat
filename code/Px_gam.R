Px.gam <- function(x,shape,rate) {
	scale <- 1/rate
	n <- 2*shape
	b0 <- scale*shape
	
	t1 <- gamma((n+1)/2)/gamma(n/2)
	t2 <- sqrt(b0/(pi*n))
	t3 <- (1 + (b0*x^2)/n)^-((n+1)/2)
	
	t1*t2*t3
}

# a little recurssion!!!!
PPx.gam <- function(x,shape,rate,comp=TRUE) {
	if(length(x)==1) {
		intgral <- integrate(Px.gam,0,x,shape=shape,rate=rate)
		if(intgral$message != "OK") print(intrgral$message)
		val <- intgral$value
		if(comp) {
			return(1 - 2*val)
		} else {
			return(2*val)
		}
	} else {
		return(sapply(x,function(X) PPx.gam(X,shape,rate,comp)))
	}
}


