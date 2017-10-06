library(paleobioDB)

## function to take a taxon and break it down into subtaxa (recursively) if
## the parent is deemed too large
#' @param x is the output from `pbdb_taxa`
#' @param maxSize is the maximum allowable clade size

subTaxa <- function(x, maxSize = 400) {
    ## clean up and exclude vertebrates
    x <- x[x$status == 'belongs to', ]
    if('Vertebrata' %in% x$taxon_name) {
        x <- x[x$taxon_name != 'Vertebrata', ]
    }
    
    ## clades to further decompose 
    tooBig <- which(x$size > maxSize)
    
    ## get children of large clades
    chil <- lapply(as.character(x$taxon_no[tooBig]), function(tname) {
        print(tname)
        out <- try(pbdb_taxa(id = tname, vocab = 'pbdb', show = 'size', 
                             rel = 'children'))
        
        if(class(out) == 'try-error') {
            if(grepl('port 80' %in% attr(out, 'condition'))) {
                Sys.sleep(10)
                out <- pbdb_taxa(id = tname, vocab = 'pbdb', show = 'size', 
                                 rel = 'children')
            } else {
                stop(attr(out, 'condition'))
            }
        }
        
        ## hack for a bug in paleobioDB where empty columns are not returned, see:
        ## https://github.com/ropensci/paleobioDB/issues/18
        if(!all(names(x) %in% names(out))) {
            out[, names(x)[!(names(x) %in% names(out))]] <- NA
        }
        out <- out[, names(x)]
        
        return(out)
    })
    
    res <- rbind(x[-tooBig, ], do.call(rbind, chil))
    
    ## recursively break-up clades
    if(any(res$size > maxSize)) {
        subTaxa(res, maxSize = maxSize)
    } else {
        return(res)
    }
}

## function to loop over taxa and retrieve their occurence data
#' @param taxaDF is the `data.frame` returned by `subTaxa`
#' @param show are all the desired additional columns

getOccs <- function(taxaDF, show) {
    ## get all the potential column names
    url <- sprintf('http://paleobiodb.org/data1.1/occs/list.txt?taxon_name=Canis&show=%s&limit=1', 
                   paste0(show, collapse = ','))
    rawAPI <- names(read.csv(url))
    
    ## loop over taxa, getting occurrences data
    dat <- lapply(1:nrow(taxaDF), function(i) {
        temp <- try(pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[i],
                                     vocab = 'pbdb', show = show))
        
        ## deal with possible errors
        if('try-error' %in% class(temp)) {
            if(grepl('C stack', attr(temp, 'condition'))) {
                ## too big
                newTaxa <- subTaxa(allTaxa[i, ], maxSize = round(allTaxa$size[i] / 2))
                temp <- getOccs(newTaxa, show = show)
            } else if(grepl('reg_count != df_count', attr(temp, 'condition'))) {
                ## nothing there, usually the `id` works
                temp <- pbdb_occurrences(limit = 'all', id = allTaxa$taxon_no[i],
                                         vocab = 'pbdb', 
                                         show = c('ident', 'phylo', 'lith', 'loc', 'time', 
                                                  'geo', 'stratext'))
            } else if(grepl('port 80' %in% attr(temp, 'condition'))) {
                ## too busy, wait a (10) sec
                Sys.sleep(10)
                temp <- pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[i],
                                         vocab = 'pbdb', 
                                         show = c('ident', 'phylo', 'lith', 'loc', 'time', 
                                                  'geo', 'stratext'))
            } else {
                stop(attr(temp, 'condition'))
            }
        }
        
        ## make sure all columns are present
        if(!any(rawAPI %in% names(temp))) {
            temp[, rawAPI[!(rawAPI %in% names(temp))]] <- NA
        }
        temp <- temp[, rawAPI]
        
        return(temp)
    })
    
    return(do.call(rbind, dat))
}

## get all taxa of interest
allTaxa <- pbdb_taxa(name = 'Animalia', vocab = 'pbdb', show = 'size', rel = 'children')
allTaxa <- subTaxa(allTaxa, maxSize = 400)
allTaxa$taxon_name <- as.character(allTaxa$taxon_name)
allTaxa <- allTaxa[order(allTaxa$size, decreasing = TRUE), ]










## loop over queries, getting occurrences data
dat <- pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[1], vocab = 'pbdb', 
                        show = show)

length(rawAPI)
length(dat)

lapply(2:nrow(allTaxa), function(i) {
    temp <- try(pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[i],
                                 vocab = 'pbdb', 
                                 show = c('ident', 'phylo', 'lith', 'loc', 'time', 
                                          'geo', 'stratext')))
    if('try-error' %in% class(temp)) {
        if(grepl('C stack', attr(temp, 'condition'))) {
            newTaxa <- subTaxa(allTaxa[, i])
        } else if(grepl('reg_count != df_count', attr(temp, 'condition'))) {
            temp <- pbdb_occurrences(limit = 'all', id = allTaxa$taxon_no[i],
                                     vocab = 'pbdb', 
                                     show = c('ident', 'phylo', 'lith', 'loc', 'time', 
                                              'geo', 'stratext'))
        } else if(grepl('port 80' %in% attr(temp, 'condition'))) {
            Sys.sleep(10)
            temp <- pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[i],
                                     vocab = 'pbdb', 
                                     show = c('ident', 'phylo', 'lith', 'loc', 'time', 
                                              'geo', 'stratext'))
        } else {
            stop(attr(temp, 'condition'))
        }
    }
})


foo <- pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[1], vocab = 'pbdb', 
                        show = c('phylo', 'lith', 'lithext', 'loc', 'time'))


