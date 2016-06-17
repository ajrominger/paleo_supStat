########  methods for sstat objects  ########
source("~/R_functions/logPlot.R")
source("~/R_functions/logAxis.R")

##	gives the likelihood under sstat model
sstat.lik <- function(par,dat) {
	this.shape <- par[1]
	this.rate <- par[2]
	
	-sum(log(Px.gam(dat,this.shape,this.rate)))
}
	
##	find the maximum likelihood estimate for sstat parameters
mle.sstat <- function(x,use.all=TRUE) {
	if(use.all) {
		these.dat <- unlist(x$Px.raw)
	} else {
		these.dat <- unlist(x$Px.sub)
	}
	
	optim(c(0.55,0.17),sstat.lik,method="BFGS",hessian=TRUE,dat=these.dat)
}


##	bootstrap likelihood
boot.mle.sstat <- function(x,B=1000,use.all=TRUE) {
	if(use.all) {
		these.dat <- x$Px.raw
	} else {
		these.dat <- x$Px.sub
	}
	
	boots <- replicate(B,
				{
					sub.dat <- sapply(these.dat,sample,size=1)
					this.mle <- try(optim(c(0.55,0.17),sstat.lik,method="BFGS",hessian=TRUE,dat=sub.dat),silent=TRUE)
					if(class(this.mle) != "try-error") {
						if(this.mle$convergence != 0) {
							out <- rep(NA,2)
						} else {
							out <- this.mle$par
						}
					} else {
						out <- rep(NA,2)
					}
					out
				})
	sstat.out <- rbind(quantile(boots[1,],c(0.025,0.975),na.rm=TRUE),quantile(boots[2,],c(0.025,0.975),na.rm=TRUE))
	rownames(sstat.out) <- c("shape","rate")
	
	norm.boots <- replicate(B, {
		sub.dat <- sapply(these.dat,sample,size=1)
		
		sd(sub.dat)
	})
	norm.out <- rbind(0,quantile(norm.boots,c(0.025,0.975),na.rm=TRUE))
	rownames(norm.out) <- c("mean","sd")
	
	return(list(sstat=sstat.out, norm=norm.out))
}


##	logLik for sstat object
logLik.sstat <- function(x,fitted=FALSE,use.all=TRUE) {
	if(use.all) {
		these.dat <- unlist(x$Px.raw)
	} else {
		these.dat <- unlist(x$Px.sub)
	}
	
	lik <- sum(log(x$Px(these.dat)))
	
	if(fitted) {
		attr(lik,"df") <- 2
	} else {
		attr(lik,"df") <- 0
	}
	
	class(lik) <- "logLik"
	
	return(lik)
}

##	AIC for sstat object (maybe defunct now?)
AIC.sstat <- function(x,use.all=TRUE) {
	if(use.all) {
		these.dat <- unlist(x$Px.raw)
	} else {
		these.dat <- unlist(x$Px.sub)
	}
	
	lik <- x$Px(these.dat)
	
	return(-sum(log(lik)))
}

##	plot for sstat object
plot.sstat <- function(x,sstat.col="red",norm.col="blue",show.norm=TRUE,add.legend=TRUE, xaxt='s',yaxt='s', y.lim, ...) {
	ecdf.all <- my.ecdf(abs(unlist(x$Px.raw)),complement=TRUE)
	ecdf.sub <- my.ecdf(abs(unlist(x$Px.sub)),complement=TRUE)
	if(missing("y.lim")) y.lim <- c(ecdf.all[,2],ecdf.sub[,2])
	
	plot(ecdf.all,log="xy",col="gray",pch=16,
		 # min of all could be 0, so make smart x-limits to avoid log(0)
		 xlim=c(ifelse(min(ecdf.all[,1])==0,min(ecdf.sub[,1]),min(ecdf.all[,1],ecdf.sub[,1])),max(ecdf.all[,1],ecdf.sub[,1])),
		 ylim=range(y.lim), xaxt='n', yaxt='n',
		 ...)
	points(ecdf.sub, cex=ifelse('cex' %in% names(list(...)), list(...)$cex, 1))
	
	if(xaxt != 'n') logAxis(1)
	if(yaxt != 'n') logAxis(2)
	
	PPx <- x$PPx
	curve(PPx(x,comp=TRUE),col=sstat.col,lwd=2,add=TRUE)
	
	if(show.norm) {
		this.sd <- sd(unlist(x$Px.sub))
		curve(2*pnorm(x,0,this.sd,lower.tail=FALSE),col=norm.col,lwd=2,add=TRUE)
	}
	
	if(add.legend) {
		leg <- c("All data",paste("Data (n >= ",x$minN,")",sep=""),"Superstatistics")
		col <- c("gray", par('fg'), sstat.col)
		pch <- c(16,1,1)
		# pt.lwd <- c(0,0,0)
		# pt.cex <- c(0,0,0)
		# lwd <- c(2,2,1)
		pt.lwd <- c(1,1,NA)
		pt.cex <- c(1,1,NA)
		lwd <- c(NA,NA,2)
		
		if("panel.first" %in% names(list(...))) {
			leg <- c(leg,"Superstatistics likelihood CI")
			col <- c(col,hsv(s=0.5))
			pt.lwd <- c(pt.lwd,1)
			pt.cex <- c(pt.cex,2)
			lwd <- c(lwd,NA)
			pch <- c(pch,15)
		}
		
		if(show.norm) {
			leg <- c(leg,"ML normal")
			col <- c(col,norm.col)
			pt.lwd <- c(pt.lwd,NA)
			pt.cex <- c(pt.cex,NA)
			lwd <- c(lwd,2)
			pch <- c(pch,1)
			
			if("panel.first" %in% names(list(...))) {
				leg <- c(leg,"Normal likelihood CI")
				col <- c(col,hsv(0.6, s=0.5))
				pt.lwd <- c(pt.lwd,1)
				pt.cex <- c(pt.cex,2)
				lwd <- c(lwd,NA)
				pch <- c(pch,15)
			}
		}
		
		extracex <- ifelse('cex' %in% names(list(...)), list(...)$cex, 1)
		legend("bottomleft",legend=leg,col=col,pch=pch,pt.lwd=pt.lwd,pt.cex=pt.cex*extracex,lwd=lwd,bty="n")
	}
}

##	function to add confidence interval polygon
##	from ML analysis
mle.poly <- function(ci.ma,fun,from=10^-2,to=10^2,...) {
	these.x <- 10^c(seq(par("usr")[1],par("usr")[2],length=25),seq(par("usr")[2],par("usr")[1],length=25))
	these.y <- c(fun(these.x[1:25],ci.ma[1,1],ci.ma[2,2]),fun(these.x[26:50],ci.ma[1,2],ci.ma[2,1]))
	these.y[these.y < 10^par("usr")[3]] <- 10^par("usr")[3]
	
	polygon(x=these.x,y=these.y,...)
}

