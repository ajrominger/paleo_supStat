logPlot <- function(x,y,log="x",...) {
	if(missing(y)) {
		if(is.null(ncol(x))) {
			y <- x
			x <- 1:length(y)
		} else {
			y <- x[,2]
			x <- x[,1]
		}
	}
	plot(x,y,log=log,axes=FALSE,...)
	usr <- par("usr")
	
	other.par <- list(...)
	yes.x <- ifelse(length(other.par$xaxt != "n") > 0, other.par$xaxt != "n", TRUE)
	yes.y <- ifelse(length(other.par$yaxt != "n") > 0, other.par$yaxt != "n", TRUE)
	
	if(yes.x) {
		if(length(grep("x",log)) > 0) {
			maj.minX <- ceiling(usr[1])
			maj.maxX <- floor(usr[2])
			
			axis(1,xaxp=c(10^maj.minX,10^maj.maxX,1))
			
			if(maj.minX - usr[1] >= log(2,10)) {
				min.incX <- rep(log(2:9,10),length(maj.minX:maj.maxX)+1)
				min.incX <- min.incX + unlist(lapply((maj.minX-1):maj.maxX,rep,times=8))
				min.incX <- min.incX[min.incX >= usr[1] & min.incX <= usr[2]]
			} else {
				min.incX <- rep(log(2:9,10),length(maj.minX:maj.maxX))
				min.incX <- min.incX + unlist(lapply(maj.minX:maj.maxX,rep,times=8))
				min.incX <- min.incX[min.incX >= usr[1] & min.incX <= usr[2]]
			}
			
			axis(1,at=10^min.incX,labels=FALSE,col.ticks="black",tcl=-0.3)
			
		} else {
			axis(1)
		}
	}
	
	if(yes.y) {
		if(length(grep("y",log)) > 0) {
			maj.minY <- ceiling(usr[3])
			maj.maxY <- floor(usr[4])
			
			axis(2,yaxp=c(10^maj.minY,10^maj.maxY,1))
			
			if(maj.minY - usr[3] >= log(2,10)) {
				min.incY <- rep(log(2:9,10),length(maj.minY:maj.maxY)+1)
				min.incY <- min.incY + unlist(lapply((maj.minY-1):maj.maxY,rep,times=8))
			} else {
				min.incY <- rep(log(2:9,10),length(maj.minY:maj.maxY))
				min.incY <- min.incY + unlist(lapply(maj.minY:maj.maxY,rep,times=8))
			}
			
			axis(2,at=10^min.incY,labels=FALSE,col.ticks="black",tcl=-0.3)
			
		} else {
			axis(2)
		}
	}
	
	box()
}
