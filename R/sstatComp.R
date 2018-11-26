sstatComp <- function(grp.data,minN=15,xlab="Absolute Fluctuation",
                      ylab="Cumulative Density",leg=TRUE,plotit=TRUE) {
	these2use <- sapply(grp.data,length) >= minN
	p2use <- grp.data[these2use]
	
	cat("computing Gaussian fit for p_k(x|sigma) \n")
	pk.par <- sapply(p2use,function(x) unlist(normLS(x)[c("par","value")]))
	pk.par <- t(pk.par)
	colnames(pk.par) <- c("mu","sig","ss")
	
	cat("re-centering \n")
	for(i in 1:length(p2use)) {
		p2use[[i]] <- p2use[[i]] - pk.par[i,"mu"]
	}
	
	cat("computing f(beta) \n")
	f.beta.par <- gammaLS(1/(pk.par[,"sig"])^2)$par
	fuent.par <- c(n=2*f.beta.par[1],b0=f.beta.par[1]*f.beta.par[2])
	
	cat("computing P(x) \n")
	this.Px <- function(x) Px.gam(x,f.beta.par[1],f.beta.par[2])
	this.PPx <- function(x,comp=TRUE) PPxGam(x,f.beta.par[1],f.beta.par[2],comp)
	
	out <- list(gam.par=f.beta.par,sspar=fuent.par,beta=1/(pk.par[,"sig"])^2,
	            sumSq=pk.par[,"ss"],minN=minN,raw.pk=p2use,
				Px.raw=grp.data,Px.sub=p2use,incld=these2use,Px=this.Px,PPx=this.PPx)
	
	class(out) <- "sstat"
	
	if(plotit) plot(out)
	
	return(out)
}

