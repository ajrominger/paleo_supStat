library(paleobioDB)

setwd('~/Dropbox/Research/paleo_supStat/data')

show <- c('ident', 'phylo', 'lith', 'loc', 'time', 'geo', 'stratext')
uri <- sprintf(
    'https://paleobiodb.org/data1.2/occs/list.csv?base_name=Animalia&show=%s&limit=all', 
    paste0(show, collapse = ',')
)

## get pbdb occurences
x <- read.csv(uri, as.is = TRUE)

## remove crinates
crinates <- pbdb_taxa(name = 'Craniata', vocab = 'pbdb', rel = 'children')
for(i in 1:nrow(crinates)) {
    x <- x[x[[crinates$rank[i]]] != crinates$taxon_name[i], ]
}

## only well lithified specimens
x <- x[x$lithification1 == 'lithified', ]

## only marine environments
x <- x[x$environment %in% readLines('marine_env.txt'), ]
