## after loading `sepkoski_data_obs.RData'

sepkoski$phy <- sapply(strsplit(as.character(sepkoski$phy)," "),function(x) x[1])
unique(sepkoski$phy)

class.info <- t(sapply(strsplit(levels(sepkoski$class)," "),function(x)
					{
						if(length(x) > 1) {
							return(c(x[1],paste(x[2:length(x)],collapse=" ")))
						} else {
							return(c(x,""))
						}
					}))