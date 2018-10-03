#' @description plot method for sstat class
#' @param x the `sstat` object
#' @param sstatCol color for super stats fit
#' @param normCol color for Gaussian fit
#' @param showNorm logical, should Gaussian fit be shown
#' @param addLegend logical, should legend be added
#' @param ... other parameters passed to `plot.default`

plot.sstat <- function(x, sstatCol = 'red', normCol = 'blue', 
                       showNorm = TRUE, addLegend = TRUE, ...) {
    thisECDF <- my.ecdf(abs(unlist(x$Px.sub)), complement = TRUE)
    
    # helper function to deal with optional axis arguments
    .axissetup <- function(side) {
        if(sprintf('%saxt', side) %in% names(pargs)) {
            if(pargs[[sprintf('%saxt', side)]] == 'n') {
                assign(sprintf('%saxfun', side), function(...) {}, pos = 1)
            } else {
                if(side %in% pargs$log) {
                    assign(sprintf('%saxfun', side), socorro::logAxis, pos = 1)
                } else {
                    assign(sprintf('%saxfun', side), axis, pos = 1)
                }
            }
        } else {
            if(side %in% pargs$log) {
                assign(sprintf('%saxfun', side), socorro::logAxis, pos = 1)
            } else {
                assign(sprintf('%saxfun', side), axis, pos = 1)
            }
        }
    }
    
    pargs <- list(...)
    if(!('log' %in% names(pargs))) pargs$log <- 'xy'
    .axissetup('x')
    .axissetup('y')
    
    pargs$xaxt <- 'n'
    pargs$yaxt <- 'n'
    
    do.call(plot, c(list(x = thisECDF), pargs))
    
    PPx <- x$PPx
    curve(PPx(x, comp = TRUE), col = sstatCol, lwd = 2, add = TRUE)
    
    if(showNorm) {
        thisSD <- sd(unlist(x$Px.sub))
        curve(2*pnorm(x, 0, thisSD, lower.tail = FALSE), col = normCol, lwd = 2, add = TRUE)
    }
    
    if(addLegend) {
        leg <- c('Observed', 'Superstatistics')
        col <- c(par('fg'), sstatCol)
        pch <- c(ifelse('pch' %in% names(list(...)), list(...)$pch, 1), NA)
        pt.lwd <- c(1, NA)
        pt.cex <- c(1, NA)
        lwd <- c(NA, 2)
        
        if('panel.first' %in% names(list(...))) {
            leg <- c(leg, 'Superstatistics likelihood CI')
            col <- c(col, socorro::colAlpha(sstatCol, 0.5))
            pt.lwd <- c(pt.lwd, 1)
            pt.cex <- c(pt.cex, 2)
            lwd <- c(lwd, NA)
            pch <- c(pch, 15)
        }
        
        if(showNorm) {
            leg <- c(leg, 'ML Gaussian')
            col <- c(col, normCol)
            pt.lwd <- c(pt.lwd, NA)
            pt.cex <- c(pt.cex, NA)
            lwd <- c(lwd, 2)
            pch <- c(pch, NA)
        }
        
        extracex <- ifelse('cex' %in% names(list(...)), list(...)$cex, 1)
        legend('bottomleft', legend = leg, col = col, pch = pch, pt.lwd = pt.lwd, 
               pt.cex = pt.cex*extracex, lwd = lwd, bty = 'n')
    }
}


#' @description function to add confidence interval polygon from ML analysis
#' @param 

mle.poly <- function(ci.ma,fun,from=10^-2,to=10^2,...) {
    these.x <- 10^c(seq(par("usr")[1],par("usr")[2],length=25),seq(par("usr")[2],par("usr")[1],length=25))
    these.y <- c(fun(these.x[1:25],ci.ma[1,1],ci.ma[2,2]),fun(these.x[26:50],ci.ma[1,2],ci.ma[2,1]))
    these.y[these.y < 10^par("usr")[3]] <- 10^par("usr")[3]
    
    polygon(x=these.x,y=these.y,...)
}
