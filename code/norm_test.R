source("~/R_functions/my_ecdf.R")

mult.shapiro.test <- function(x, plotit=TRUE, ...) {
	stats <- sapply(x, function(X) shapiro.test(X)$p.value)
	
	if(plotit) {
		plot(my.ecdf(stats), xlim = 0:1, ylim = 0:1, xlab="P-values", ylab="Cumulative density", ...)
		abline(0,1,col="red")
	}
	browser()
	ks.test(stats, punif, alternative = "g")
}

##  error analysis
# tru.false <- replicate(1000, {
	# true <- replicate(400,shapiro.test(rnorm(10+rpois(1,20),0,rgamma(1,0.65,0.8)))$p.value)
	# false <- replicate(400,shapiro.test(rlnorm(10+rpois(1,20),10,0.2))$p.value)
	
	# c(ks.test(true, punif, alternative = "g")$p.value, ks.test(false, punif, alternative = "g")$p.value)
# })

# sum(tru.false[2,] < 0.05)
