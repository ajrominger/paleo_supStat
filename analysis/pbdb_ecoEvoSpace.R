library(socorro)
library(plyr)

pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)

head(pbdbDat[, 30:37])

foo <- ddply(pbdbDat, 'class', function(x) {
    data.frame(v = sum(!duplicated(x[, c(30, 32:37)])))
})

hist(foo$v[-1], breaks = seq(0, max(foo$v[-1]) + 1, by = 1))

names(pbdbDat)
??'nonmetric'
