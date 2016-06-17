##	needed data
load("/Users/andrewrominger/Desktop/Research/paleo_supStat/data/sepk/sepkoski_GOOD.RData")
sepk.date <- rev(read.table("/Users/andrewrominger/Desktop/Research/paleo_supStat/data/sepk/sepkoski_invert.txt",header=TRUE)$date)[-1]


##	needed functions

# convenience function for making data.frames of fluctuation
stage.fact <- as.factor(colnames(sepkoski[,-(1:6)]))
make.flux.tab <- function(X) {
	lapply(X,function(x)
		{
			x <- x[,-(1:2)]
			div <- apply(x,2,sum)
			flux <- diff(c(0,div))
			fl <- apply(x,1,function(xx)
						{
							these.occ <- which(xx==1)
							c(min(these.occ),max(these.occ))
						})
			frst <- table(stage.fact[fl[1,]])[stage.fact]
			last <- table(stage.fact[fl[2,]])[stage.fact]
			data.frame(date=sepk.date,div=div,flux=flux,frst=frst,last=last)
		})
}

# function for computing ``events'' i.e. ignors time steps w/o change
event.flux <- function(x) {
	events <- which(diff(x) != 0)
	diff(x[c(events[1],events+1)])
}

##	make data.frame containing only certain orders

#sepk.uncert <- sepkoski[grep("ORD_UNCERT",sepkoski$ord),]
sepk.cert <- sepkoski[-grep("ORD_UNCERT",sepkoski$ord),]
sepk.cert <- sepk.cert[-which(sepk.cert$ord %in% c("?CONULARIIDA","CERATITIDA","CHELONIELLIDA")),]


