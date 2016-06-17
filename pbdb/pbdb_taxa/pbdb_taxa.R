these.ord <- names(pbdb.sstat.ord.cor$beta)
these.cls <- as.character(pbdb.dat$occurrences.class_name[match(these.ord,pbdb.dat$occurrences.order_name)])


pbdb.taxa <- read.csv("~/Desktop/Research/paleo_supStat/pbdb/pbdb_taxa/opinions.csv",header=TRUE,stringsAsFactors=FALSE)

pbdb.rank <- read.csv("~/Desktop/Research/paleo_supStat/pbdb/pbdb_taxa/valid_taxa.csv",header=TRUE,stringsAsFactors=FALSE)
pbdb.rank <- pbdb.rank[pbdb.rank$taxon_rank %in% c("class","superorder",
												   "subclass","phylum","infraorder",
												   "infraclass","superclass","subphylum",
												   "superkingdom","kingdom","subkingdom",
												   "superphylum"),]
rownames(pbdb.rank) <- NULL
pbdb.rank[(1:3)+nrow(pbdb.rank),1:7] <- rbind(c("class","Demospongia","A. J.","Rominger","","",""),
											  c("subclass","Protobranchia","A. J.","Rominger","","",""),
											  c("class","Echinus","A. J.","Rominger","","",""))

pbdb.rank[nrow(pbdb.rank)-(0:2),8] <- 2012


pbdb.taxa <- pbdb.taxa[pbdb.taxa$child_name %in% pbdb.rank$taxon_name | pbdb.taxa$parent_name %in% pbdb.rank$taxon_name,]

pbdb.hier <- pbdb.taxa[pbdb.taxa$status == "belongs to",]
pbdb.synm <- pbdb.taxa[pbdb.taxa$status != "belongs to",]

##	get rid of multiple parents
taxa.web <- with(pbdb.hier,tapply(parent_name,child_name, function(x) {
	if(length(x) > 1) {
		# check if any are misspelled, etc
		if(any(x %in% pbdb.synm$child_name)) {
			x <- sapply(x,function(X) {
				if(X %in% pbdb.synm$child_name) {
					return(pbdb.synm$parent_name[pbdb.synm$child_name==X])
				} else {
					return(X)
				}
			
			})
		}
		
		# get rid of duplicates from misspelling, etc.
		x <- unique(x)
		
		# if still multiple parents
		if(length(x) > 1) {
			# take most recent opinion
			x <- x[order(pbdb.hier$pubyr[match(x,pbdb.hier$child_name)])][1]
		}
	}
	
	return(x)
}))

taxa.tree <- cbind(unlist(taxa.web),names(taxa.web))
rownames(taxa.tree) <- NULL

##	get rid of "animalia"
taxa.tree <- taxa.tree[taxa.tree[,1] != "Animalia",]

##	recursive function to find all parents
find.parent <- function(x,rank.default="order") {
	taxa.ord <- c("kingdom","subkingdom","superphylum",
				  "phylum","subphylum","superclass",
				  "class","subclass","infraclass",
				  "superorder","order")
	
	name2add <- pbdb.rank$taxon_rank[pbdb.rank$taxon_name==x[1]]
	if(length(name2add) == 0) {
		if(length(x) == 1) {
			name2add <- rank.default
		} else {
			name2add <- "unknown"
		}
	}
	names(x) <- c(name2add,names(x)[-1])
	
	parent <- taxa.tree[taxa.tree[,2]==x[1],1]
	old.names <- names(x)
	
	if(length(parent) > 0) {
		x <- c(parent,x)
		if(x[1] == x[2]) {
			out <- x[taxa.ord]
			names(out) <- taxa.ord
			return(out)
		} else {
			return(find.parent(x))
		}
	} else {
		out <- x[taxa.ord]
		names(out) <- taxa.ord
		return(out)
	}
}


taxon.hier <- lapply(these.ord[these.ord %in% taxa.tree],find.parent)
taxon.hier <- do.call(rbind,taxon.hier)

##	a few fixes
taxon.hier <- taxon.hier[,!apply(taxon.hier,2,function(x) all(is.na(x)))]
taxon.hier <- taxon.hier[!is.na(taxon.hier[,"order"]),]
taxon.hier <- taxon.hier[,colnames(taxon.hier) != "superclass"]

##	missing taxa
taxon.hier[taxon.hier[,"order"] == "Neoloricata",c("kingdom","phylum")] <- c("Metazoa","Mollusca")
taxon.hier[taxon.hier[,"order"] == "Camarodonta",c("kingdom","phylum","class","superorder")] <- c("Metazoa","Echinodermata","Echinoidae","Echinacea")

all(taxon.hier[,"order"] %in% names(pbdb.sstat.ord.cor$beta))

these.cls[!(names(pbdb.sstat.ord.cor$beta) %in% taxon.hier[,"order"])] %in% taxon.hier[,"class"]


this.phy <- taxon.hier[match(names(pbdb.sstat.ord.cor$beta),taxon.hier[,"order"]),"phylum"]
phy.var.obs <- mean(tapply(log(pbdb.sstat.ord.cor$beta),this.phy,var))

phy.var.sim <- replicate(1000,{
	mean(tapply(sample(log(pbdb.sstat.ord.cor$beta),rep=FALSE),this.phy,var))
})

plot(density(phy.var.sim))
abline(v=phy.var.obs)
abline(v=quantile(phy.var.sim,c(0.025,0.975)),lty=2)






##	make a phylogeny
library(ade4)
library(ape)

pbdb.taxo <- as.data.frame(taxon.hier[,ncol(taxon.hier):1])
rownames(pbdb.taxo) <- pbdb.taxo[,1]
pbdb.taxo <- pbdb.taxo[,-1]

as.taxo(pbdb.taxo[,-7])












!(taxa.tree[,2] %in% pbdb.rank$taxon_name[pbdb.rank$taxon_rank %in% 
											c("order","class","superorder",
												   "subclass","phylum","infraorder",
												   "infraclass","superclass","subphylum",
												   "superkingdom","kingdom","subkingdom",
												   "superphylum")])


all(sapply(taxa.web,length) == 1)