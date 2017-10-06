library(paleobioDB)

subTaxa <- function(x) {
    ## max taxon size allowable
    maxSize <- 400
    
    ## clean up and exclude vertebrates
    x <- x[x$status == 'belongs to', ]
    if('Vertebrata' %in% x$taxon_name) {
        x <- x[x$taxon_name != 'Vertebrata', ]
    }
    
    
    ## clades to further decompose 
    tooBig <- which(x$size > maxSize)
    
    ## get children of large clades
    chil <- lapply(as.character(x$taxon_name[tooBig]), function(tname) {
        print(tname)
        out <- try(pbdb_taxa(name = tname, vocab = 'pbdb', show = 'size', 
                             rel = 'children'))
        
        if(class(out) == 'try-error') {
            if(grepl('port 80' %in% attr(out, 'condition'))) {
                Sys.sleep(10)
                out <- pbdb_taxa(name = tname, vocab = 'pbdb', show = 'size', 
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
        
        ## avoid circular taxonomy, an issue with, e.g. Nautilida and Nautiloidea
        ## being each other's parent
        if(any(x$parent_no %in% out$taxon_no)) {
            ## take out the problem child
            out <- out[!(out$taxon_no %in% x$parent_no), ]
            
            ## make sure we don't keep looking
            x$size[x$parent_no %in% out$taxon_no] <<- 1
            
            ## make sure not to remove the parent from `x`
            tooBig <<- tooBig[as.character(x$taxon_name[tooBig]) != tname]
        }
        
        return(out)
    })
    
    res <- rbind(x[-tooBig, ], do.call(rbind, chil))
    
    ## recursively break-up clades
    if(any(res$size > maxSize)) {
        subTaxa(res)
    } else {
        return(res)
    }
}

allTaxa <- pbdb_taxa(name = 'Animalia', vocab = 'pbdb', show = 'size', rel = 'children')
allTaxa <- subTaxa(allTaxa)
allTaxa$taxon_name <- as.character(allTaxa$taxon_name)

foo <- pbdb_occurrences(limit = 'all', base_name = allTaxa$taxon_name[191], vocab = 'pbdb', 
                        show = c('phylo', 'lith', 'lithext', 'loc', 'time'))
head(foo)

