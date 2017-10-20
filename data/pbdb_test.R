library(paleobioDB)
library(parallel)

setwd('~/Dropbox/Research/paleo_supStat/data')

## call to the API
show <- paste0(c('ident', 'phylo', 'lith', 'loc', 'time', 'geo', 'stratext', 
                 'ecospace'), 
               collapse = ',')
version <- '1.2'
base_name <- 'Neopronorites'
min_ma <- 2.5
max_ma <- 550
timerule <- 'contain'
envtype <- 'marine'

uri <- sprintf(
    'https://paleobiodb.org/data%s/occs/list.csv?base_name=%s&show=%s&limit=all&min_ma=%s&max_ma=%s&timerule=%s&envtype=%s', 
    version,
    base_name,
    show,
    min_ma,
    max_ma,
    timerule, 
    envtype
)

y <- read.csv(uri, as.is = TRUE)

write.csv(y, file = 'test.csv', row.names = FALSE)
