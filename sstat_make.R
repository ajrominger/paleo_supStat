# **script to put everything together**

# get PBDB data through API and clean
source('data/pbdb_data_get.R')

# make a taxonomic hash
source('data/pbdb_taxa.R')

# correct for incomplete sampling and bias
source('data/pbdb_3TPub_make.R')

# compare diversity curves
source('analysis/pbdb_divCurve.R')

# conduct super-statistical analysis
source('analysis/pbdb_sstat.R')

# permute taxa and analyze null fit of super-stats to data
source('analysis/pbdb_dperm.R')

# explore occupation of eco-space by taxa
source('analysis/pbdb_ecoEvoSpace.R')

# explore relationship between beta and richness
source('analysis/pbdb_betaRichness.R')
