#' @description plot method for sstat class
#' @param x the `sstat` object
#' @param sstatCol color for super stats fit
#' @param normCol color for Gaussian fit
#' @param showNorm logical, should Gaussian fit be shown
#' @param addLegend logical, should legend be added
#' @param ... other parameters passed to `plot.default`

plot.sstat <- function(x, sstatCol = 'red', normCol = 'blue', 
                       showNorm = TRUE, addLegend = TRUE, ...) {
    thisECDF <- socorro::simpECDF(abs(unlist(x$Px.sub)), complement = TRUE)
    
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
    if(!('xlab' %in% names(pargs))) pargs$xlab <- '|Fluctuations|'
    if(!('ylab' %in% names(pargs))) pargs$ylab <- 'Cumulative density'
    .axissetup('x')
    .axissetup('y')
    
    pargs$xaxt <- 'n'
    pargs$yaxt <- 'n'
    
    do.call(plot, c(list(x = thisECDF), pargs))
    xaxfun(1)
    yaxfun(2)
    
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
            leg <- c(leg, 'Superstatistics CI')
            col <- c(col, socorro::colAlpha(sstatCol, 0.25))
            pt.lwd <- c(pt.lwd, 1)
            pt.cex <- c(pt.cex, 2)
            lwd <- c(lwd, NA)
            pch <- c(pch, 15)
        }
        
        if(showNorm) {
            leg <- c(leg, 'Gaussian')
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
#' @param ci the matrix of CI intervals for the parameter values returned by `bootMLE.sstat`
#' @param fun the CDF function to plug the parameter values into
#' @param ... further arguments passed to `polygon` (e.g. `col`, `boarder`, etc.)

mlePoly <- function(ci, fun, ...) {
    n <- 50
    x <- seq(par('usr')[1], par('usr')[2], length = n)
    x <- c(x, rev(x))
    
    if(par('xlog')) x <- 10^x
    
    y <- c(fun(x[1:n], ci[1, 1], ci[2, 2]), fun(x[(1:n) + n], ci[1, 2], ci[2, 1]))
    
    polygon(x = x, y = y, ...)
}
