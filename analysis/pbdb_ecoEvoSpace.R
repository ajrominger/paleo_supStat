# **a script to evaluate hour occupancy of eco-evolutionary space changes 
# across taxonomy**

pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)

# extract only the eco/evo/life history data and remove duplicates
# `taxon_environment`, `reproduction`, `ontogeny`
eeSpaceDat <- pbdbDat[, c('phylum', 'class', 'order', 'family', 'otu', 
                          'taxon_environment', 'motility', 'life_habit',
                          'vision', 'diet', 'reproduction', 'ontogeny')]
eeSpaceDat <- eeSpaceDat[!duplicated(eeSpaceDat), ]

# remove entries that are all missing
eeSpaceDat <- eeSpaceDat[rowSums(eeSpaceDat[, -(1:5)] != '') != 0, ]


#' function to determine how many eco-evo hypercubes are occupied by each taxonomic level
#' @param tax the taxonomic unit to consider
#' @param eeDat a data.frame containing eco-evo data
eeOcc <- function(tax, eeDat) {
    sapply(split(eeDat[tax != '', ], tax[tax != '']), 
           function(x) sum(!duplicated(x)))
}

#' bootstraps ee occupancy to compare with permutations
#' @param x the vector of niche occupancies
#' @param B number of bootrap replicates
#' @param fun the function to apply to each replicate
eeOccBoot <- function(x, B, fun) {
    replicate(B, fun(sample(x, length(x), replace = TRUE)))
}


#' function to permute genera and calculate null occupancy
#' @param tax the taxonomic unit to consider
#' @param eeDat a data.frame containing eco-evo data
#' @param B number of bootrap replicates
#' @param fun the function to apply to each replicate
eeOccPerm <- function(tax, eeDat, B, fun) {
    replicate(B, {
        eeDat <- eeDat[sample(nrow(eeDat), size = nrow(eeDat), replace = TRUE), ]
        fun(eeOcc(tax, eeDat))
    })
}


famEE <- eeOccBoot(eeOcc(eeSpaceDat$family, eeSpaceDat[, -(1:5)]), 500, mean)
ordEE <- eeOccBoot(eeOcc(eeSpaceDat$order, eeSpaceDat[, -(1:5)]), 500, mean)
clsEE <- eeOccBoot(eeOcc(eeSpaceDat$class, eeSpaceDat[, -(1:5)]), 500, mean)
phyEE <- eeOccBoot(eeOcc(eeSpaceDat$phylum, eeSpaceDat[, -(1:5)]), 500, mean)

famEEPerm <- eeOccPerm(eeSpaceDat$family, eeSpaceDat[, -(1:5)], 500, mean)
ordEEPerm <- eeOccPerm(eeSpaceDat$order, eeSpaceDat[, -(1:5)], 500, mean)
clsEEPerm <- eeOccPerm(eeSpaceDat$class, eeSpaceDat[, -(1:5)], 500, mean)
phyEEPerm <- eeOccPerm(eeSpaceDat$phylum, eeSpaceDat[, -(1:5)], 500, mean)



plot(1:4, ylim = range(famEE, phyEEPerm), type = 'n', log = 'y')
segments(x0 = 1:4, y0 = c(min(famEE), min(ordEE), min(clsEE), min(phyEE)), 
         y1 = c(max(famEE), max(ordEE), max(clsEE), max(phyEE)))
segments(x0 = 1:4, y0 = c(min(famEEPerm), min(ordEEPerm), min(clsEEPerm), min(phyEEPerm)), 
         y1 = c(max(famEEPerm), max(ordEEPerm), max(clsEEPerm), max(phyEEPerm)))

points(1:4, c(mean(famEE), mean(ordEE), mean(clsEE), mean(phyEE)), pch = 16)
points(1:4, c(mean(famEEPerm), mean(ordEEPerm), mean(clsEEPerm), mean(phyEEPerm)), pch = 16)





hist(foo)






# extract unique niches
eeUnique <- unique(eeSpaceDat[, -(1:5)])

# map unique niches to all rows of `eeSpaceDat`
eeUniqueName <- apply(eeUnique, 1, paste, collapse = ';')
eeSpaceDat$nicheID <- match(apply(eeSpaceDat[, -(1:5)], 1, paste, collapse = ';'), 
                            eeUniqueName)

# map each unique niche to all the others that could subsume it (due to missing data)
eeShort <- gsub('-', '', eeUniqueName)
sapply(gsub('-', '.*', eeUniqueName), function(x) {
    out <- rep(0, nrow(eeUnique))
    out[grep(x, eeShort)] <- 1
    return(out)
})


# loop through unique niches and classify which are truely unique and which should be
# subsummed due to missing data
eeUnique <- unique(foo)
eeUnique <- eeUnique[order(rowSums(is.na(eeUnique)), decreasing = TRUE), ]

eeMap <- matrix(0, nrow = nrow(eeUnique), ncol = nrow(eeUnique))
diag(eeMap) <- 1


for(i in 1:max(which(rowSums(is.na(eeUnique)) > 0))) {
    x <- unlist(eeUnique[i, ])
    theseNA <- is.na(x)
    x <- x[!theseNA]
    dat <- eeUnique[1:nrow(eeUnique) > i, !theseNA, drop = FALSE]
    
    # these life histories subsume life history i
    subsume <- apply(dat, 1, function(r) all(x == r))
    eeMap[i, which(subsume)] <- 1
}











sapply(1:max(which(rowSums(is.na(eeUnique)) > 0)), function(i) {
    x <- unlist(eeUnique[i, ])
    theseNA <- is.na(x)
    x <- x[!theseNA]
    dat <- eeUnique[1:nrow(eeUnique) > i, !theseNA]
}) 





# make life histories numers
eeSpaceDat[, -(1:5)] <- apply(eeSpaceDat[, -(1:5)], 2, function(x) as.integer(as.factor(x)))

head(eeSpaceDat)
sum(eeSpaceDat$vision == '')

foo <- rbind(c(NA, 'a', 'b', 'c'), 
             c('1', 'a', 'b', 'c'), 
             c(NA, 'x', 'b', NA), 
             c(NA, 'x', 'b', NA), 
             c(NA, 'z', 'b', 'c'), 
             c('2', 'z', NA, 'c'))
# foo[is.na(foo)] <- '-'

bla1 <- apply(foo, 1, paste, collapse = ';')
bla2 <- gsub('-', '.*', bla1)
sapply(bla2, function(x) grepl(x, bla1))


regexpr('abc', 'a.*c')


removeLooseDupRows <- function(x) { 
    # browser()
    if (nrow(x) <= 1) 
        return(x) 
    ii <- do.call("order", 
                  args=lapply(seq_len(ncol(x)), 
                              function(col) x[ , col])) 
    
    dup_index <- logical(nrow(x)) 
    i0 <- -1 
    for (k in 1:length(ii)) { 
        i <- ii[k] 
        if (any(is.na(x[i, ]))) { 
            if (i0 == -1) 
                next 
            if (any(x[i, ] != x[i0, ], na.rm=TRUE)) 
                next 
            dup_index[i] <- TRUE 
        } else { 
            i0 <- i 
        } 
    } 
    
    dup_index
}

foo
foo[!removeLooseDupRows(foo) | duplicated(foo), ]




foo <- foo[order(rowSums(is.na(foo))), 
           order(colSums(is.na(foo)), decreasing = TRUE)]


cbind(foo, cbind(duplicated(foo), 
      sapply(1:(ncol(foo) - 1), function(i) duplicated(foo[, -(1:i)]))))


apply(foo, 2, duplicated)

calcSpace <- function(x) {
    browser()
    
    rmNACol <- function(m) {
        m[, colSums(is.na(m)) == 0]
    }
    
    rowSums(is.na(x))
    
    
    
    x <- x[order(rowSums(!is.na(x)), decreasing = TRUE), ]
    dups <- sapply(1:(ncol(x) - 1), function(i) duplicated(x[, -(1:i)]))
    dups <- cbind(duplicated(x), dups)
    dups %*% t(1 * (!is.na(x)))
}
calcSpace(foo)



foo <- ddply(pbdbDat, 'otu', function(x) {
    data.frame(v = sum(!duplicated(x[, c(30, 32:37)])))
})

hist(foo$v[-1], breaks = seq(0, max(foo$v[-1]) + 1, by = 1))

names(pbdbDat)
??'nonmetric'
