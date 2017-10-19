options(stringsAsFactors = FALSE)

setwd('~/Dropbox/Research/paleo_supStat/data')

coreURI <- 'http://fossilworks.org/bridge.pl?a=displayInterval&interval_no='

tbinInfo <- lapply(1:1108, function(i) {
    print(i)
    linfo <- try(readLines(paste0(coreURI, i), n = 150))
    
    if('try-error' %in% class(linfo)) 
        linfo <- try(readLines(paste0(coreURI, i), n = 150))
    
    if('try-error' %in% class(linfo)) {
        thisTbin <- thisMax <- thisMin <- NA
    } else {
        thisTbin <- gsub('^.*10 million year bin: |<br>.*$', '', 
                         linfo[grep('10 million year bin', linfo)])
        thisMax <- as.numeric(gsub('^.*Lower boundary: equal to | Ma.*$|[^0-9\\.]', '',
                                   linfo[grep('Lower boundary: equal to', linfo)]))
        thisMin <- as.numeric(gsub('^.*Upper boundary: equal to | Ma.*$|[^0-9\\.]', '',
                                   linfo[grep('Upper boundary: equal to', linfo)]))
    }
    
    return(data.frame(tbin = ifelse(length(thisTbin) == 0, NA, thisTbin), 
                      ma_min = ifelse(length(thisMin) == 0, NA, thisMin),
                      ma_max = ifelse(length(thisMax) == 0, NA, thisMax)))
})

tbinInfo <- do.call(rbind, tbinInfo)

tbins <- plyr::ddply(tbinInfo, 'tbin', function(x) range(x[, -1], na.rm = TRUE))
tbins <- tbins[!is.na(tbins$tbin), ]
names(tbins) <- c('tbin', 'ma_min', 'ma_max')
tbins <- tbins[order(tbins$ma_min), ]

write.csv(tbins, 'tbins2.csv', row.names = FALSE)
