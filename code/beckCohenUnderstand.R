library(MASS)

k <- 20
X <- replicate(1000,sum(rnorm(k,sd=1/sqrt(rgamma(1,0.6,0.8)))^2))
X <- replicate(1000,sum(rnorm(k,sd=0.2)^2))

plot(ecdf(X),xlim=c(min(X),1000),log="x")
curve(pchisq(x,df=k-10),col="red",add=TRUE)
curve(pgamma(x, shape=k/2, scale=2*mean(X)/k),col="blue",add=TRUE)

mean(X)

##  SO...
##  chisq = gamma when var of norms = 1. when vars != 1 (but still
##  var_i = var_j \forall i,j) then sum of norms^2 is still
##  gamma. this seems related to fact that scale = 2/k for case where
##  vars = 1 and when vars != 1 we need a scale, specifically <x>
##  which is not accounted for by the chisq. when var_i != var_j then
##  everything breaks down

par <- fitdistr(X,dgamma,list(shape=1,scale=1))
with(par,curve(pgamma(x,estimate[1],estimate[2]),col="blue",lty=2,add=TRUE))

X <- replicate(1000,sum(rnorm(k,sd=1/sqrt(rgamma(1,0.6,0.8)))^2))
plot(ecdf(k/X),log="x",xlim=c(min(k/X),10))

X <- replicate(1000,var(rnorm(k,sd=1/sqrt(rgamma(1,0.6,0.8)))))
plot(ecdf(1/X),log="x",xlim=c(min(1/X),10))

curve(pgamma(x,0.6,0.8),col="red",lty=2,add=TRUE)


