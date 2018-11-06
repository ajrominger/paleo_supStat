# **a script to compare our 3TPub curve to other estimates of richness through the Phanerozoic**

# package with diversity dynamics subsampling functions
library(divDyn)

# load and prep data
pbdbDat <- read.csv('data/pbdb_data.csv', as.is = TRUE)
tbin <- read.csv('data/tbinsMid.csv', as.is = TRUE)
tbin$tbin <- factor(tbin$tbin, levels = tbin$tbin)
pbdbDat$tbin <- factor(pbdbDat$tbin, levels = levels(tbin$tbin))
pbdbDat$tbinNum <- as.integer(pbdbDat$tbin)

pbdbDatUnique <- pbdbDat[!duplicated(paste0(pbdbDat$collection_no, pbdbDat$otu)), ]



data(corals)
data(stages)
fossils <- corals[corals$stg != 95, ]

# indicate identical collection/genus combinations
collGenus <- paste(fossils$collection_no, fossils$genus)
# omit the duplicates from the occurrence datasets
fossGen <- fossils[!duplicated(collGenus),]


# rarifaction considering collections as well as occurrences
subUWunit <- subsample(fossGen, bin="stg", tax="genus",
                       iter=100, q=10, type="cr", unit="collection_no")
# classic rarifaction
subCR <- subsample(fossGen, bin="stg", tax="genus", iter=100, q=40)

# occurrence-weighted by list
subOW <- subsample(fossGen, bin="stg", tax="genus", coll="collection_no",
                   iter=100, q=40, type="oxw")



# sqs with good's u based on different singletons
sqsCollSing <-subsample(fossGen, iter=50, q=0.4,
                        tax="genus", bin="stg", type="sqs",
                        singleton="occ")
sqsRefSing <-subsample(fossGen, iter=50, q=0.4,
                       tax="genus", bin="stg", type="sqs", ref="reference_no",
                       singleton="ref")
sqsNoSing <-subsample(fossGen, iter=50, q=0.4,
                      tax="genus", bin="stg", type="sqs", singleton=FALSE)

# sqs with collection correction
sqsByColl <-subsample(fossGen, iter=50, q=0.5,
                      tax="genus", bin="stg", ref="reference_no",coll="collection_no",
                      type="sqs", singleton="ref", excludeDominant=T, largestColl=T, byList=TRUE)

tsplot(stages, shading="series", boxes="per", xlim=51:95,
       ylab="corrected SIB diversity", ylim=c(-3, 3.5))
lines(stages$mid[1:94], scale(subCR$divCSIB), col="black")
# lines(stages$mid[1:94], scale(subUWunit$divCSIB), col="blue")
# lines(stages$mid[1:94], scale(subOW$divCSIB), col="red")
lines(stages$mid[1:94], scale(sqsByColl$divCSIB), col="blue")
