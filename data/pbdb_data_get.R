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

## remove unnecceary columns
c2rm <- c('record_type', 'reid_no', 'flags', 'identified_name', 
          'identified_rank', 'identified_no', 'difference', 'species_name', 
          'species_reso', 'lithdescript', 'lithology1', 'minor_lithology1', 
          'lithology2', 'lithification2', 'minor_lithology2', 'cc', 'state', 
          'county', 'latlng_basis', 'geogcomments', 'geology_comments', 
          'zone', 'localsection', 'localbed', 'localorder', 
          'regionalsection', 'regionalbed', 'regionalorder', 
          'stratcomments')
x <- x[, !(names(x) %in% c2rm)]

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


## standard time bins
tbins <- read.csv('timebins.csv', as.is = TRUE)
x$tbin <- ''
for(i in 1:nrow(tbins)) {
    x$tbin[x$max_ma <= tbins$ma_max[i] & x$min_ma >= tbins$ma_min[i]] <- tbins$tbin[i]
}

x <- x[x$tbin != '', ]

## write out fully processed data
write.csv(x, 'pbdb_data.csv', row.names = FALSE)
