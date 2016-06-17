source("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdb_sqs.R")

ammonite <- read.csv("/Users/andrewrominger/Desktop/Research/paleo_supStat/pbdbTest-occs.csv")

amm.sqs.par <- calc.sqs.par(ammonite)
names(amm.sqs.par)

amm.sqs <- pbdb.sqs(ammonite,q=0.5,pars=amm.sqs.par,B=2)


bla <- calc.sqs(ammonite,0.8,amm.sqs.par,raw.id=TRUE)

bla <- pbdb.sqs(ammonite,0.8,amm.sqs.par,B=2,raw.id=TRUE)

bla.lev <- unique(as.character(ammonite$occurrences.family_name[ammonite$collection_no %in% unique(unlist(bla$ID))]))


sapply(bla$ID,function(x)
						{
							x <- factor(x,levels=bla.lev)
							
						})



##	maybe best way is to create a vector of orders that can be used with colXtaxa
##	the taxa present will then be associated with an order
##	one problem is how to easily track orders across time bins...if taxa are factors
##	then all taxa will always be present in colXtaxa...that seems to be the case!!!

dim(calc.sqs.par(ammonite[1:100,])$colXtaxa)

length(unique(ammonite$collection_no))