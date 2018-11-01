#' helper function to calculate corrected flux
#' @param x the matrix of corrected diversities over which to calculate fluxes

calcFlux <- function(x) {
    apply(x, 2, function(X) {
        flux <- diff(c(0, X))
        return(flux[flux != 0])
    })
}
