# **script to run super stat analysis on PBDB data**

# source needed functions
R.utils::sourceDirectory('R', modifiedOnly = FALSE)
library(socorro)

# load and prepare data
# ---------------------

pbdbFamDiv <- read.csv('data/pbdb_3TPub-corrected.csv', row.names = 1)

# coarsen to higher taxonomic groupings


# tbin midpoints
tbinMid <- read.csv('data/tbinsMid.csv', as.is = TRUE)
tbinNames <- tbinMid$tbin
tbinMid <- as.numeric(tbinMid[, 2])
names(tbinMid) <- tbinNames


# super stat analysis for families
# --------------------------------

# corrected flux
pbdbFamFlux <- apply(pbdbFamDiv, 2, function(x) {
    flux <- diff(c(0, x))
    return(flux[flux != 0])
})

# make sstat object
sstatPBDBord3TP <- sstatComp(pbdbFamFlux, minN = 10, plotit = FALSE)

# likelihood CI for family-level sstat analysis
sstatPBDBord3TPCI <- bootMLE.sstat(sstatPBDBord3TP, B = 1000, useAll = FALSE)

# plot it
pdf('ms/fig_Px.pdf', width = 4.25, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.5, 0), cex.lab = 1.4)
plot(sstatPBDBord3TP, xlim = c(1e-04, 4e+01), ylim = c(8e-05, 1),
     xaxt = 'n', yaxt = 'n',
     panel.first = quote(mlePoly(sstatPBDBord3TPCI$sstat, PPx.gam,
                                 col = hsv(alpha = 0.25), border = NA)))
logAxis(1:2, expLab = TRUE)

dev.off()


# plot p_k(x|b) and f(beta) for families
# --------------------------------------

# idea for normality test: sample 1 from each order and do ks test on that subsampled set

# highlight these families

# make CDF for all scale family-level fluctuations
pAll <- lapply(sstatPBDBord3TP$raw.pk, 
               function(x) simpECDF(scale(x)[, 1], complement = TRUE))
pAll <- do.call(rbind, pAll)


layout(matrix(c(1, 2, 1, 3), nrow = 2))

par(oma = c(0, 3, 0, 0) + 0.25, mar = c(4, 0, 0, 0) + 0.25, 
    mgp = c(2, 0.5, 0), cex.lab = 1.4)


# low = Tainoceratidae

pbdbFamDiv[, 'Tainoceratidae']


# mid = Lophospiridae
# hi = Chonetidae

plot(1, xlim = c(540, 0), xaxt = 'n', xaxs = 'i', xlab = '')
paleoAxis(1)
mtext('Millions of years ago', side = 1, line = 3.5)

par(mar = c(3, 0, 1, 0) + 0.25)
plot(pAll, xlim = c(-4, 4), col = 'gray', ylim = c(0, 1.025),
     xlab = 'Scaled fluctuations', ylab = 'Cumulative density')
curve(pnorm(x, lower.tail = FALSE), lwd = 2, add = TRUE)

plot(simpECDF(sstatPBDBord3TP$beta, complement = TRUE), ylim = c(0, 1.025),
     log = 'x', xaxt = 'n', yaxt = 'n',
     xlab = expression(beta), col = 'gray')
logAxis(1, expLab = TRUE)
curve(pgamma(x, sstatPBDBord3TP$gam.par[1], sstatPBDBord3TP$gam.par[2], 
             lower.tail = FALSE), 
      col = 'black', lwd = 2, add = TRUE)

# low = Tainoceratidae
# mid = Lophospiridae
# hi = Chonetidae

low <- which(sstatPBDBord3TP$beta < 2^-(2 - 0.25) & sstatPBDBord3TP$beta > 2^-(2 + 0.25))
mid <- which(sstatPBDBord3TP$beta < 2^(0.5 + 0.25) & sstatPBDBord3TP$beta > 2^(0.5 - 0.25))
hi <- which(sstatPBDBord3TP$beta < 2^(3 + 0.25) & sstatPBDBord3TP$beta > 2^(3 - 0.25))

sort(colSums(pbdbFamDiv[, names(hi)] > 0))


