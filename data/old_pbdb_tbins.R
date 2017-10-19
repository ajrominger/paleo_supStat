foo <- read.csv('pbdb_2013-05-28/marInv-occs.csv', as.is = TRUE)


bla <- plyr::ddply(foo[foo$collections.10_my_bin != '', ], 
                   'collections.10_my_bin', function(dat) {
    c(min = min(dat$ma_min), max = max(dat$ma_max))
})

names(bla) <- c('tbin', 'ma_min', 'ma_max')

bla <- bla[order(bla$ma_min), ]

for(i in 2:nrow(bla)) {
    bla[i, 2] <- bla[i-1, 3] <- mean(c(bla[i, 2], bla[i-1, 3]))
}
bla[1, 2] <- 0

write.csv(bla, 'timebins.csv', row.names = FALSE)
