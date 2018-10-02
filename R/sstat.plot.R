#	plot for sstat object
plot.sstat <- function(x, sstat.col = 'red', norm.col = 'blue', show.norm = TRUE, add.legend = TRUE, xaxt = 's', yaxt = 's', y.lim, ...) {
    ecdf.all <- my.ecdf(abs(unlist(x$Px.raw)), complement = TRUE)
    ecdf.sub <- my.ecdf(abs(unlist(x$Px.sub)), complement = TRUE)
    if(missing('y.lim')) y.lim <- c(ecdf.sub[, 2], ecdf.sub[, 2])
    
    plot(ecdf.sub, log = 'xy', col = 'gray', pch = 16, 
         # min of all could be 0, so make smart x-limits to avoid log(0)
         xlim = c(ifelse(min(ecdf.all[, 1]) == 0, min(ecdf.sub[, 1]), min(ecdf.all[, 1], ecdf.sub[, 1])), max(ecdf.all[, 1], ecdf.sub[, 1])), 
         ylim = range(y.lim), xaxt = 'n', yaxt = 'n', 
         ...)
    # points(ecdf.sub, cex = ifelse('cex' %in% names(list(...)), list(...)$cex, 1))
    
    if(xaxt !=  'n') logAxis(1)
    if(yaxt !=  'n') logAxis(2)
    
    PPx <- x$PPx
    curve(PPx(x, comp = TRUE), col = sstat.col, lwd = 2, add = TRUE)
    
    if(show.norm) {
        this.sd <- sd(unlist(x$Px.sub))
        curve(2*pnorm(x, 0, this.sd, lower.tail = FALSE), col = norm.col, lwd = 2, add = TRUE)
    }
    
    if(add.legend) {
        leg <- c('All data', paste('Data (n >=  ', x$minN, ')', sep = ''), 'Superstatistics')
        col <- c('gray', par('fg'), sstat.col)
        pch <- c(16, 1, 1)
        # pt.lwd <- c(0, 0, 0)
        # pt.cex <- c(0, 0, 0)
        # lwd <- c(2, 2, 1)
        pt.lwd <- c(1, 1, NA)
        pt.cex <- c(1, 1, NA)
        lwd <- c(NA, NA, 2)
        
        if('panel.first' %in% names(list(...))) {
            leg <- c(leg, 'Superstatistics likelihood CI')
            col <- c(col, hsv(s = 0.5))
            pt.lwd <- c(pt.lwd, 1)
            pt.cex <- c(pt.cex, 2)
            lwd <- c(lwd, NA)
            pch <- c(pch, 15)
        }
        
        if(show.norm) {
            leg <- c(leg, 'ML normal')
            col <- c(col, norm.col)
            pt.lwd <- c(pt.lwd, NA)
            pt.cex <- c(pt.cex, NA)
            lwd <- c(lwd, 2)
            pch <- c(pch, 1)
            
            if('panel.first' %in% names(list(...))) {
                leg <- c(leg, 'Normal likelihood CI')
                col <- c(col, hsv(0.6, s = 0.5))
                pt.lwd <- c(pt.lwd, 1)
                pt.cex <- c(pt.cex, 2)
                lwd <- c(lwd, NA)
                pch <- c(pch, 15)
            }
        }
        
        extracex <- ifelse('cex' %in% names(list(...)), list(...)$cex, 1)
        legend('bottomleft', legend = leg, col = col, pch = pch, pt.lwd = pt.lwd, pt.cex = pt.cex*extracex, lwd = lwd, bty = 'n')
    }
}