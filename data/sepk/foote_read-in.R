######  Basic read-in and manipulation  ######
this.dir <- "/Users/andrewrominger/Desktop/Fulbright/P1_Macroevolution/Data/Sepkoski"
file <- paste(this.dir,"compendium_pri_mod.txt",sep="/")
foote.check <- readLines(file)
foote.check <- foote.check[-(1:max(grep("#",foote.check)))]   # get rid of meta-data
foote.check <- foote.check[foote.check != ""]	# get rid of empty rows

######  get rid of synonyms  #######
foote.check <- foote.check[-grep("@",foote.check)]
foote.check <- foote.check[-which(substr(foote.check,1,1)=="=")]

######  manage higher taxonomy  ######
require(zoo)

nlines <- length(foote.check)
sepk.phy <- rep(NA,nlines)		# make 'empty' vectors for storage
sepk.cls <- rep(NA,nlines)
sepk.ord <- rep(NA,nlines)

# Phyla
these.phy <- grep("Ph. ",foote.check)
sepk.phy[these.phy] <- gsub("Ph. ","",foote.check[these.phy])
sepk.phy <- na.locf(sepk.phy,na.rm=FALSE,fromLast=FALSE)
sepk.phy <- sapply(strsplit(sepk.phy," "),function(x) x[1])

# Classes
these.cls <- grep("Cl. ",foote.check)
# note: sometimes no class or order assignment, use this to catch...
sepk.cls[these.phy] <- paste(sepk.phy[these.phy],"CLS_UNCERT",sep="_")
sepk.cls[these.cls] <- gsub("Cl. ","",foote.check[these.cls])
sepk.cls <- na.locf(sepk.cls,na.rm=FALSE,fromLast=FALSE)
sepk.cls <- sapply(strsplit(sepk.cls," "),function(x) x[1])



# Orders
these.ord <- grep("Or. ",foote.check)
sepk.ord[these.phy] <- paste(sepk.phy[these.phy],"ORD_UNCERT",sep="_")
sepk.ord[these.cls] <- paste(sepk.cls[these.cls],"ORD_UNCERT",sep="_")
sepk.ord[these.ord] <- gsub("Or. ","",foote.check[these.ord])
sepk.ord <- na.locf(sepk.ord,na.rm=FALSE,fromLast=FALSE)
sepk.ord <- sapply(strsplit(sepk.ord," "),function(x) x[1])

# figure out what's what (print and search online)
#unique(sepk.phy)
#unique(sepk.cls[sepk.phy=="NEMATODA"])		# bilaterian, (check for terrestrial)
#unique(sepk.cls[sepk.phy=="ANNELIDA"])		# bilaterian, (check for terrestrial)
#unique(sepk.cls[sepk.phy=="PROBLEMATICA"])	# incertae sedis
#unique(sepk.cls[sepk.phy=="MOLLUSCA"])		# bilaterian, (check for terrestrial)
#unique(sepk.cls[sepk.phy=="LOBOPODA"])		# super group containing ARTHROPODA?!?!
#unique(sepk.cls[sepk.phy=="ARTHROPODA"])	# bilaterian, (check for terrestrial)
#unique(sepk.ord[sepk.phy=="LOBOPODA"])
#unique(sepk.ord[sepk.cls=="ARACHNIDA"])
#unique(sepk.ord[sepk.cls=="BRANCHIOPODA"])
#unique(sepk.ord[sepk.cls=="MALACOSTRACA"])
#unique(sepk.ord[sepk.cls=="OSTRACODA"])
#unique(sepk.ord[sepk.cls=="GASTROPODA"])
#unique(sepk.ord[sepk.cls=="BIVALVIA"])
# all other checking needs doing on genera....

# now clean-up all vectos
sepk.tax <- cbind(sepk.phy,sepk.cls,sepk.ord)[-c(these.phy,these.cls,these.ord),]
foote.check <- foote.check[-c(these.phy,these.cls,these.ord)]

######  mange genera  #######
# get rid of ++ (but store in vector if we want them)
nopls <- numeric(length(foote.check))
nopls[substr(foote.check,1,1)==" "] <- 1
substr(foote.check,1,1) <- ""	# gets rid of ++

# take first 'word' and that's genus!
sepk.gen <- substr(foote.check,1,regexpr(" ",foote.check)-1)

# check genera....NOTE: bad organization, but must be done before trimming
# foote.check above (under `now clean-up all vectors')
#length(sepk.gen) == length(sepk.phy)
#unique(sepk.gen[sepk.ord=="ISOPODA"])
#unique(sepk.gen[sepk.ord=="PODOCOPIDA"])
#unique(sepk.gen[sepk.ord=="BASOMMATOPHORA"])
#unique(sepk.gen[sepk.cls=="BRANCHIOPODA"])
#unique(sepk.gen[sepk.phy=="LOBOPODA"])
#unique(sepk.gen[sepk.phy=="NEMATODA"])
#unique(sepk.gen[sepk.cls=="COPEPODA"])

nogo.phy <- c("CHORDATA","ACTINOPODA","CILIOPHORA","RHIZOPODEA","PORIFERA","TRILOBOZOA")
nogo.cls <- c("COPEPODA","UNIONOIDA","MYRIAPODA")


######  mange stratigraphy  ######
##	convenience function to take matrix
##	of first/last and produce 'life table'
make.life.tab <- function(X,mar=1) {
	res <- apply(X,mar,function(xx)
			{
				xx[xx==0] <- NA
				x1 <- na.locf(xx,na.rm=FALSE,fromLast=FALSE)
				x2 <- na.locf(xx,na.rm=FALSE,fromLast=TRUE)
				xx <- x1 + x2
				xx[is.na(xx)] <- 0
				xx/2
			})
	t(res)
}

# extract first/last stage (if they exist)
sepk.stage <- read.table("~/Desktop/Fulbright/P1_Macroevolution/Data/Sepkoski/sepkoski_invert.txt",header=TRUE)$stage
sepk.stage <- sepk.stage[-c(1,length(sepk.stage))] # take off pre-Cm and Recent (for now)
sepk.stage <- as.character(sepk.stage)			# make character
sepk.stage <- paste("\\(",sepk.stage,sep="")	# for better matching
sepk.stage <- c("- R",sepk.stage)				# add back in Recent (in form suited for matching)
sepk.stage <- rev(sepk.stage)					# time order

##	object holding first/last occurances
sepk.read.fl <- sapply(sepk.stage,function(x)
						{
							1*(regexpr(x,foote.check) > 0)
						#	num.match <- gregexpr(x,foote.check)
						#	sapply(num.match,function(X) sum(X>0))
						})

# which are poorely resolved--update all objects accordingly
bad.strat <- apply(sepk.read.fl,1,sum) < 1
sepk.read.fl <- sepk.read.fl[!bad.strat,]
sepk.tax <- sepk.tax[!bad.strat,]
sepk.gen <- sepk.gen[!bad.strat]
nopls <- nopls[!bad.strat]
foote.check <- foote.check[!bad.strat]

# evaluate uncertanties...i.e. '?'
uncert.strat <- rep(1,length(sepk.gen))
uncert.strat[grep("\\?",foote.check)] <- 0	# 0 for uncertain

# interpolate between first and last...and clean up colnames...
substr(colnames(sepk.read.fl),1,2) <- ""
substr(colnames(sepk.read.fl),1,1) <- ""	# note sure why need to do 2x...

sepk.read.fl[sepk.read.fl==2] <- 1	# make everything '1'
sepk.read.fl <- make.life.tab(sepk.read.fl)

sepkoski <- data.frame(sepk.tax,nopls,sepk.gen,uncert.strat,sepk.read.fl)
colnames(sepkoski)[c(1:3,5)] <- c("phy","class","ord","gen")
for(i in 1:3) {
	sepkoski[,i] <- as.character(sepkoski[,i])
}

sepkoski.all <- sepkoski
sepkoski <- sepkoski[-which(sepkoski$phy %in% nogo.phy),]
sepkoski <- sepkoski[-which(sepkoski$class %in% nogo.cls),]

###### SAVE IT !!!!!
save(sepkoski,sepkoski.all,file="~/Desktop/sepkoski_GOOD.RData")

## make objects to hold phy/clade/ord flux

# convenience function for making data.frames of fluctuation
stage.fact <- as.factor(colnames(sepkoski[,-(1:6)]))
make.flux.tab <- function(X) {
	lapply(X,function(x)
		{
			x <- x[,-(1:2)]
			div <- apply(x,2,sum)
			flux <- c(0,diff(div))
			fl <- apply(x,1,function(xx)
						{
							these.occ <- which(xx==1)
							c(min(these.occ),max(these.occ))
						})
			frst <- table(stage.fact[fl[1,]])[stage.fact]
			last <- table(stage.fact[fl[2,]])[stage.fact]
			data.frame(date=sepk.date,div=div,flux=flux,frst=frst,last=last)
		})
}

sepk.phy <- split(sepkoski[,5:ncol(sepkoski)],sepkoski$phy)
sepk.cls <- split(sepkoski[,5:ncol(sepkoski)],sepkoski$class)
sepk.ord <- split(sepkoski[,5:ncol(sepkoski)],sepkoski$ord)


sepk.phy.divflux <- make.flux.tab(sepk.phy)
sepk.cls.divflux <- make.flux.tab(sepk.cls)
sepk.ord.divflux <- make.flux.tab(sepk.ord)

### save these....
fi <- "/Users/andrewrominger/Desktop/sepk_superstat"
save(sepk.phy.divflux,sepk.cls.divflux,sepk.ord.divflux,file=paste(fi,"sepk_sstat_obs.RData",sep="/"))


##	what's up with Orders AMMONOIDEA, CERATITIDA and CHEILOSTOMATA?
sepkoski[sepkoski$ord=="AMMONOIDEA",1:5]	# sub-class (1324 genera)
sepkoski[sepkoski$ord=="CERATITIDA",1:5]	# ``order'' within AMMONOIDEA (426 genera)
sepkoski[sepkoski$ord=="CHEILOSTOMATA",1:5]	# kitchen sink? (575 genera)





######### all the old checking still there.....

#####  some checking
raw.foote <- readLines(file)
all(levels(sepkoski$phy) %in% gsub("Ph. ","",raw.foote[grep("Ph. ",raw.foote)]))
all(gsub("Ph. ","",raw.foote[grep("Ph. ",raw.foote)]) %in% levels(sepkoski$phy))
gsub("Ph. ","",raw.foote[grep("Ph. ",raw.foote)])[!(gsub("Ph. ","",raw.foote[grep("Ph. ",raw.foote)]) %in% levels(sepkoski$phy))]

all(levels(sepkoski$class) %in% gsub("Cl. ","",raw.foote[grep("Cl. ",raw.foote)]))
all(gsub("Cl. ","",raw.foote[grep("Cl. ",raw.foote)]) %in% levels(sepkoski$class))
gsub("Cl. ","",raw.foote[grep("Cl. ",raw.foote)])[!(gsub("Cl. ","",raw.foote[grep("Cl. ",raw.foote)]) %in% levels(sepkoski$class))]

all(levels(sepkoski$ord) %in% gsub("Or. ","",raw.foote[grep("Or. ",raw.foote)]))
all(gsub("Or. ","",raw.foote[grep("Or. ",raw.foote)]) %in% levels(sepkoski$ord))

all(gsub("Or. ","",raw.foote[grep("Or. ",raw.foote)]) %in% unique(sepk.ord[!is.na(sepk.ord)]))

gsub("Or. ","",raw.foote[grep("Or. ",raw.foote)])[!(gsub("Or. ","",raw.foote[grep("Or. ",raw.foote)]) %in% levels(sepkoski$ord))]


###### PLOT IT !!!!!
source("~/R_functions/paleoPlot.R")
sepk.div <- apply(sepkoski[,-(1:6)],2,sum)
sepk.div.con <- apply(sepkoski[,-(1:6)]*sepkoski$uncert.strat,2,sum)
sepk.date <- rev(read.table("~/Desktop/Fulbright/P1_Macroevolution/Data/Sepkoski/sepkoski_invert.txt",header=TRUE)$date)[-1]
paleoPlot(sepk.date,sepk.div,type="l",y.lim=c(185,5500),ylab="Genus Diversity",lwd=2)
lines(sepk.date,sepk.div.con)
plot(sort(sepk.div),sepk.div.con[order(sepk.div)],type="l")
abline(0,1,lty=2)


##### little checking....
plot(sepk.phy.divflux[[1]]$div,type="s")
points(cumsum(sepk.phy.divflux[[1]]$flux),type="s",col="red")
points(cumsum(sepk.phy.divflux[[1]]$frst-c(0,sepk.phy.divflux[[1]]$last[-74])),type="s",col="blue")

plot(apply(sapply(sepk.phy.divflux,function(x) x$div),1,sum),type="s")
points(sepk.div,col="red")

plot(cumsum(apply(sapply(sepk.phy.divflux,function(x) x$flux),1,sum)),type="s")
points(sepk.div-sepk.div[1],type="s",col="red")	# take off first and OK

plot(cumsum(apply(sapply(sepk.phy.divflux,function(x) x$frst),1,sum)) -
	 c(0,cumsum(apply(sapply(sepk.phy.divflux,function(x) x$last),1,sum))[-74]),type="s")
points(sepk.div,col="red")


##### plotting....
colz <- rainbow(length(sepk.ord.divflux),end=0.85,alpha=0.8)

paleoPlot(sepk.ord.divflux[[1]]$date,sepk.ord.divflux[[1]]$div,
		  type="l",ylab="Genus diversity",col=colz[1],
		  y.lim=c(0,1150),yaxt="n",frame.plot=FALSE)
for(i in 2:length(sepk.ord.divflux)) {
	lines(sepk.ord.divflux[[i]]$date,sepk.ord.divflux[[i]]$div,col=colz[i])
}
axis(2,at=seq(0,600,200))

lines(sepk.phy.divflux[[1]]$date,(sepk.div-min(sepk.div))/15 + 850,lwd=2)
axis(2,at=seq(850,(max(sepk.div)-min(sepk.div))/15 + 850,1000/15),
	 labels=seq(0,5000,by=1000))
abline(h=850)
segments(sepk.ord.divflux[[1]]$date,-5,sepk.ord.divflux[[1]]$date,-15)






