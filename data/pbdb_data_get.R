library(paleobioDB)
oldOp <- options(stringsAsFactors = FALSE) # can't set this within pbdb functions

## function to take a taxon and break it down into subtaxa (recursively) if
## the parent is deemed too large
#' @param x is the output from `pbdb_taxa`
#' @param maxSize is the maximum allowable clade size
#' @param genMin logical, should taxa be decomposed below the rank of genus 
#' (no: `genMin = TRUE`)

subTaxa <- function(x, maxSize = 400, genMin = TRUE) {
    ## clean up and exclude vertebrates
    x <- x[x$status == 'belongs to', ]
    if('Vertebrata' %in% x$taxon_name) {
        x <- x[x$taxon_name != 'Vertebrata', ]
    }
    
    ## clades to further decompose
    if(genMin) {
        tooBig <- which(x$size > maxSize & 
                            !grepl('genus|species', x$rank, ignore.case = TRUE))
    } else {
        tooBig <- which(x$size > maxSize)
    }
    
    
    if(length(tooBig) > 0) {
        ## get children of large clades
        chil <- lapply(as.character(x$taxon_no[tooBig]), function(tname) {
            print(tname)
            out <- try(pbdb_taxa(id = tname, vocab = 'pbdb', show = 'size', 
                                 rel = 'children'))
            
            if(class(out) == 'try-error') {
                if(grepl('port 80', attr(out, 'condition'))) {
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
            ## remove synonyms
            res <- res[res$status == 'belongs to', ]
            return(res)
        }
    } else {
        ## return raw `x` if no clades were too big
        return(x)
    }
}

## function to loop over taxa and retrieve their occurence data
#' @param taxaDF is the `data.frame` returned by `subTaxa`
#' @param show are all the desired additional columns
#' @param rawAPI the raw column names from the API, without missing entries removed

getOccs <- function(taxaDF, show, rawAPI) {
    ## loop over taxa, getting occurrences data
    dat <- lapply(1:nrow(taxaDF), function(i) {
        print(taxaDF[i, c('taxon_no', 'taxon_name')])
        
        temp <- try(pbdb_occurrences(limit = 'all', base_name = taxaDF$taxon_name[i],
                                     vocab = 'pbdb', show = show), silent = TRUE)
        # browser()
        ## deal with possible errors
        if('try-error' %in% class(temp)) {
            if(grepl('C stack', attr(temp, 'condition'))) {
                ## too big
                if(grepl('gen|spe', taxaDF$rank[i], ignore.case = TRUE)) {
                    ## if rank is genus, we need to go below genus, but then include the 
                    ## original genus name (to get those not described below gen) like this:
                    ## `pbdb_occurrences(limit = 'all', taxon_name = 'Nucula')`
                    newTaxa <- subTaxa(taxaDF[i, ], maxSize = ceiling(taxaDF$size[i] / 2), 
                                       genMin = FALSE)
                    temp1 <- pbdb_occurrences(limit = 'all', taxon_name = 'Nucula', 
                                              vocab = 'pbdb', show = show)
                    temp1[, rawAPI[!(rawAPI %in% names(temp1))]] <- NA
                    temp1 <- temp1[, rawAPI]
                    
                    temp2 <- getOccs(newTaxa, show = show, rawAPI = rawAPI)
                    
                    temp <- rbind(temp1, temp2)
                } else {
                    ## otherwise subtaxa as normal
                    newTaxa <- subTaxa(taxaDF[i, ], maxSize = ceiling(taxaDF$size[i] / 2))
                    temp <- getOccs(newTaxa, show = show, rawAPI = rawAPI)
                }
            } else if(grepl('reg_count != df_count', attr(temp, 'condition'))) {
                ## nothing there, usually the `id` works
                temp <- try(pbdb_occurrences(limit = 'all', id = taxaDF$taxon_no[i],
                                             vocab = 'pbdb', 
                                             show = show), 
                            silent = TRUE)
                
                ## if failed, there's no records there
                if('try-error' %in% class(temp)) {
                    temp <- as.list(rep(NA, length(rawAPI)))
                    names(temp) <- rawAPI
                    temp <- as.data.frame(temp)
                    temp$taxon_no <- taxaDF$taxon_no[i]
                    temp$taxon_name <- taxaDF$taxon_name[i]
                }
            } else if(grepl('port 80', attr(temp, 'condition'))) {
                ## too busy, wait 10 sec
                Sys.sleep(10)
                temp <- getOccs(taxaDF[i, ], show = show, rawAPI = rawAPI)
            } else {
                temp <- as.list(rep(NA, length(rawAPI)))
                names(temp) <- rawAPI
                temp <- as.data.frame(temp)
                temp$taxon_no <- taxaDF$taxon_no[i]
                temp$taxon_name <- taxaDF$taxon_name[i]
            }
        }
        
        ## make sure all columns are present
        if(any(!(rawAPI %in% names(temp)))) {
            temp[, rawAPI[!(rawAPI %in% names(temp))]] <- NA
        }
        temp <- temp[, rawAPI]
        
        return(temp)
    })
    
    return(do.call(rbind, dat))
}

## unabashed wrapper function around `getOccs` that simply minimizes the number
## of times we have to make a request from the API
getOccsWrapper <- function(taxaDF, show) {
    ## get all the potential column names
    url <- sprintf(
        'http://paleobiodb.org/data1.1/occs/list.txt?taxon_name=Canis&show=%s&limit=1', 
        paste0(show, collapse = ','))
    rawAPI <- names(read.csv(url))
    
    getOccs(taxaDF, show, rawAPI)
}

## get all taxa of interest
allTaxa <- pbdb_taxa(name = 'Animalia', vocab = 'pbdb', show = 'size', rel = 'children')
allTaxa <- subTaxa(allTaxa, maxSize = 400)
allTaxa <- allTaxa[order(allTaxa$size, decreasing = TRUE), ]


## get occurence data
show <- c('ident', 'phylo', 'lith', 'loc', 'time', 'geo', 'stratext')
allOccs <- getOccsWrapper(allTaxa, show)


## re-set options
options(oldOp)
