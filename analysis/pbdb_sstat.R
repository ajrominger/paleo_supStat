# **script to run super stat analysis on PBDB data**

# source needed functions
R.utils::sourceDirectory('R')

# load and prepare data
# ---------------------

pbdbFamDiv <- read.csv('data/pbdb_3TPub-corrected.csv', row.names = 1)

# coarsen to higher taxonomic groupings


# super stat analysis for families
# --------------------------------

# corrected flux
pbdbFamFlux <- apply(pbdbFamDiv, 2, function(x) {
    flux <- diff(c(0, x))
    return(flux[flux != 0])
})


# sstat analysis
sstatPBDBord3TP <- sstatComp(pbdbFamFlux, minN = 10, plotit = TRUE)
