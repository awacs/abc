######################################################################
#
# misc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: normalise, namesWarningFilter, calc_distance
#
######################################################################


normalise <- function(x,y){
  if(mad(y) == 0)
    return (x)
  else
    return (x/mad(y))
}


namesWarningFilter <- function(x){
    if( any( grepl( "No parameter names are given", x) ) ) invokeRestart( "muffleWarning" )
    if( any( grepl( "No summary statistics names are given", x) ) ) invokeRestart( "muffleWarning" )
}


calc_distance <- function(sumstat, target, gwt, method = "euclidean") {
    nss <- ncol(sumstat)
    weights <- sapply(1:nss, function(j) {
        m <- mad(sumstat[, j][gwt])
        if (m == 0) return(1)
        return(1 / m)
    })
    diff <- sweep(as.matrix(sumstat), 2, target)
    wdiff <- sweep(diff, 2, weights, "*")
    if (method == "euclidean") return(sqrt(rowSums(wdiff^2)))
    if (method == "manhattan") return(rowSums(abs(wdiff)))
    if (method == "wasserstein") {
        if (any(sumstat < 0) || any(target < 0))
            warning("'wasserstein' distance expects sumstat and target to be PMFs (non-negative values over ordered bins).", call. = FALSE)
        sim_cdfs <- t(apply(as.matrix(sumstat), 1, cumsum))
        obs_cdf <- cumsum(target)
        return(rowSums(abs(sweep(sim_cdfs, 2, obs_cdf, "-"))))
    }
    stop(paste("Unknown distance method:", method))
}
