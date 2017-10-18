library(paleobioDB)
library(parallel)

setwd('~/Dropbox/Research/paleo_supStat/data')

## call to the API
show <- paste0(c('ident', 'phylo', 'lith', 'loc', 'time', 'geo', 'stratext', 
                 'ecospace'), 
               collapse = ',')
version <- '1.2'
base_name <- 'Animalia^Craniata'
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

## get pbdb occurences
x <- read.csv(uri, as.is = TRUE)

## write out raw data
write.csv(x, 'pbdb_data_raw.csv', row.names = FALSE)

## clean up

## only well lithified specimens
x <- x[x$lithification1 == 'lithified', ]

## no basin-scale collections
x <- x[!(x$geogscale %in% c('', 'basin')), ]

## fine scale stratigraphy only
x <- x[!(x$stratscale %in% c('', 'group', 'supergroup')), ]


## resolve taxonomy to genus or subgenus where availible
otu <- x$genus
otu[x$subgenus_name != ''] <- ifelse(x$subgenus_reso[x$subgenus_name != ''] == '', 
                                     x$subgenus_name[x$subgenus_name != ''], 
                                     otu[x$subgenus_name != ''])
otu[x$primary_reso != ''] <- ''
x$otu <- otu
x <- x[x$otu != '', ]

 
## combine multiple records of same otu per collection
x <- mclapply(split(x, x$collection_no), mc.cores = 6, FUN = function(coll) {
    coll[match(unique(coll$otu), coll$otu), ]
})
x <- do.call(rbind, x)

# 
# ## standard time bins
# names(x)
# head(x[, c('early_interval', 'late_interval', 'max_ma', 'min_ma', '')])
# 
# names(x)

