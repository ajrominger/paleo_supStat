# pdf for P(x) with f(beat) ~ Gamma
#' @param x diversity fluctuation value
#' @param shape the shape parameter of the gamma distribution
#' @param rate the rate parameter  of the gamma distribution

Px.gam <- PxGam <- function(x, shape, rate) {
	scale <- 1 / rate
	n <- 2 * shape
	b0 <- scale * shape
	
	t1 <- gamma((n+1) / 2) / gamma(n / 2)
	t2 <- sqrt(b0 / (pi * n))
	t3 <- (1 + (b0 * x^2) / n)^-((n + 1) / 2)
	
	t1 * t2 * t3
}

# cdf for P(x) with f(beat) ~ Gamma
#' @param x diversity fluctuation value
#' @param shape the shape parameter of the gamma distribution
#' @param rate the rate parameter  of the gamma distribution
#' @param comp logical, whether to compute the complement or not (`comp = TRUE` is 
#' equivilant to `lower.tail = FALSE` for typical `p` functions [e.g. `pnorm`])

PPx.gam <- PPxGam <- function(x, shape, rate, comp=TRUE) {
	if(length(x) == 1) {
		intgral <- integrate(PxGam, 0, x, shape = shape, rate = rate)
		if(intgral$message != 'OK') print(intrgral$message)
		val <- intgral$value
		
		if(comp) {
			return(1 - 2 * val)
		} else {
			return(2 * val)
		}
		
	} else {
	    # recursive handeling for multiple `x` values
		return(sapply(x, function(X) PPxGam(X, shape, rate, comp)))
	}
}


