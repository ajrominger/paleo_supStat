pbdbRaw <- read.csv('data/pbdb_data_raw.csv', as.is = TRUE)
pbdbNew <- read.csv('data/pbdb_data.csv', as.is = TRUE)
pbdbOld <- read.csv('data/old/pbdb_2011_07-17/marInv-occs.csv', as.is = TRUE)

dim(pbdbNew)
dim(pbdbOld)

names(pbdbNew)
names(pbdbOld)

sum(unique(pbdbOld$collection_no) %in% pbdbRaw$collection_no) / length(unique(pbdbOld$collection_no))


sum(unique(pbdbOld$collection_no) %in% x$collection_no) / length(unique(pbdbOld$collection_no))

head(pbdbOld[pbdbOld$collection_no %in% pbdbRaw$collection_no & !(pbdbOld$collection_no %in% x$collection_no), ])

head(pbdbRaw[pbdbRaw$collection_no == 177, ])
