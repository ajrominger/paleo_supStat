# **script to make a taxonomic hash**

pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)
names(pbdbDat)

pbdbTax <- as.matrix(pbdbDat[, c('phylum', 'class', 'order', 'family', 'otu')])

pbdbTax <- pbdbTax[!duplicated(pbdbTax), ]

write.csv(pbdbTax, 'data/pbdb_taxa.csv', row.names = FALSE)
