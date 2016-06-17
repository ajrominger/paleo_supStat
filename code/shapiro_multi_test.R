source("~/R_functions/ks_stat.R")

x <- replicate(400,shapiro.test(rnorm(round(runif(1,10,100)),mean=rnorm(1,0,0.5),sd=1/sqrt(rgamma(1,2,1))))$statistic)

hist(x)

ks.test(x=x,y=punif)
