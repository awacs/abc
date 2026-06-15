######################################################################
#
# corrdim.R
#
# Part of the R/abc package
# Contains: corrDim4abc, is.corrdim4abc, plot.corrdim4abc, summary.corrdim4abc
#
# Grassberger-Procaccia correlation dimension (D2) of the simulated
# summary-statistic cloud, measured in the EXACT representation and
# distance metric the abc() run uses.  D2 is the effective dimension the
# cloud occupies at scale eps, and a noise-inclusive UPPER bound on the
# data's informative dimension:  informative (<= n_params) <= D2 <= k.
#
# D2 is a property of the simulated cloud ALONE -- no observed data is
# needed.  An optional 'target' (off by default) only overlays the ABC
# operating point eps_abc and the pointwise acceptance curve C_obs(eps).
#
######################################################################


## Internal: embed sumstat into the space where plain pairwise distance
## reproduces calc_distance() for the requested method.
##   euclidean / manhattan : per-column bank-MAD weighting (== calc_distance)
##   wasserstein           : per-row CDF, then L1  (== calc_distance, no weights)
.corrdim_embed <- function(S, gwt, distance) {
    if (distance == "wasserstein") {
        embed  <- t(apply(S, 1, cumsum))
        method <- "manhattan"
    } else {
        ## weights = 1 / MAD per column over the good rows, exactly as
        ## calc_distance()/normalise() compute them.
        w <- apply(S[gwt, , drop = FALSE], 2, function(col) {
            m <- mad(col)
            if (m == 0) 1 else 1 / m
        })
        embed  <- sweep(S, 2, w, "*")
        method <- distance
    }
    list(embed = embed, method = method)
}


## Internal: local slope  D2(eps) = d log C / d log eps, evaluated on the
## full eps_grid so every call returns an aligned, equal-length vector
## (lets pooled and per-replicate slopes share one grid).
## smooth = TRUE uses the derivative of a smoothing spline of log C vs
## log eps; otherwise finite differences interpolated back onto the grid.
## C is 0 (log = -Inf) at small eps; those points are fit-excluded and the
## slope there is extrapolated -- they fall in the shaded, excluded regime.
.corrdim_slope <- function(eps_grid, C, smooth) {
    le <- log(eps_grid)
    lC <- log(C)
    ok <- is.finite(le) & is.finite(lC)
    if (sum(ok) < 4) return(rep(NA_real_, length(eps_grid)))
    if (smooth) {
        fit <- smooth.spline(le[ok], lC[ok])
        as.numeric(predict(fit, le, deriv = 1)$y)
    } else {
        mid  <- (le[-1] + le[-length(le)]) / 2
        d2   <- diff(lC) / diff(le)
        good <- is.finite(d2) & is.finite(mid)
        if (sum(good) < 2) return(rep(NA_real_, length(eps_grid)))
        approx(mid[good], d2[good], xout = le, rule = 2)$y
    }
}


corrDim4abc <- function(sumstat, target = NULL, distance = "euclidean",
                        subset = NULL, tol = 0.01, nsub = 4000, nrep = 1,
                        ngrid = 250, smooth = TRUE,
                        scaling.quantiles = c(0.05, 0.40),
                        seed = NULL, ...) {
    ## 'target' is OPTIONAL and off by default. The correlation dimension is a
    ## property of the simulated cloud alone; supplying 'target' only adds the
    ## eps_abc operating-point marker + pointwise C_obs acceptance curve.

    ## checks
    ## ######
    if (!any(distance == c("euclidean", "manhattan", "wasserstein")))
        stop("'distance' must be euclidean, manhattan or wasserstein.", call. = FALSE)
    if (!is.matrix(sumstat) && !is.data.frame(sumstat))
        stop("'sumstat' has to be a matrix or data.frame.", call. = FALSE)
    if (length(scaling.quantiles) != 2 || scaling.quantiles[1] >= scaling.quantiles[2])
        stop("'scaling.quantiles' must be two increasing values in (0,1).", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)

    S       <- as.matrix(sumstat)
    numsim  <- nrow(S)
    numstat <- ncol(S)

    statnames <- if (length(colnames(S))) colnames(S) else paste0("S", 1:numstat)

    ## good rows (drop NA) and optional user subset, as in abc()/cv4abc()
    ## #################################################################
    gwt <- rep(TRUE, numsim)
    gwt[attributes(na.omit(S))$na.action] <- FALSE
    if (is.null(subset)) subset <- rep(TRUE, numsim)
    gwt    <- as.logical(gwt * subset)
    n.good <- sum(gwt)
    if (n.good < 4) stop("Fewer than 4 usable simulations.", call. = FALSE)

    ## metric-identity embedding
    ## #########################
    emb        <- .corrdim_embed(S, gwt, distance)
    embed      <- emb$embed
    distmethod <- emb$method

    ## operating point eps_abc + pointwise acceptance curve C_obs(eps)
    ## ##############################################################
    ## Uses the package's own calc_distance() and the abc() tol->quantile
    ## rule (abc.R: nacc <- ceiling(length(dist)*tol); ds <- sort(dist)[nacc]).
    eps_abc <- NA_real_
    d_obs   <- NULL
    if (!is.null(target)) {
        target <- as.numeric(target)
        if (length(target) != numstat)
            stop("'target' length must equal the number of summary statistics.", call. = FALSE)
        d_all   <- calc_distance(S, target, gwt, method = distance)
        d_obs   <- d_all[gwt]
        nacc    <- ceiling(n.good * tol)
        eps_abc <- sort(d_obs)[nacc]
    }

    ## all-pairs correlation integral on nrep random subsamples
    ## ########################################################
    m <- min(nsub, n.good)
    good.idx <- which(gwt)
    pooled <- vector("list", nrep)
    for (r in seq_len(nrep)) {
        idx        <- sample(good.idx, m)
        pooled[[r]] <- as.numeric(dist(embed[idx, , drop = FALSE], method = distmethod))
    }
    all.d <- unlist(pooled)

    ## common log-spaced eps grid spanning the pooled pairwise distances
    pos    <- all.d[all.d > 0]
    eps_grid <- exp(seq(log(min(pos)), log(max(all.d)), length.out = ngrid))

    ## C(eps) pooled + per replicate (fraction of pairs within eps)
    C       <- ecdf(all.d)(eps_grid)
    C_reps  <- t(vapply(pooled, function(d) ecdf(d)(eps_grid), numeric(ngrid)))
    C_obs   <- if (!is.null(d_obs)) ecdf(d_obs)(eps_grid) else NULL

    ## local slope D2(eps), pooled + per replicate -- all on the common grid
    eps_mid <- eps_grid
    D2      <- .corrdim_slope(eps_grid, C, smooth)
    D2_reps <- t(vapply(seq_len(nrep),
                        function(r) .corrdim_slope(eps_grid, C_reps[r, ], smooth),
                        numeric(ngrid)))

    ## scaling region + headline plateau
    ## #################################
    qs     <- quantile(all.d, scaling.quantiles, names = FALSE)
    region <- eps_mid >= qs[1] & eps_mid <= qs[2]
    if (!any(region)) region <- rep(TRUE, length(eps_mid))
    D2_reg <- D2[region]
    plateau <- list(
        value  = median(D2_reg, na.rm = TRUE),
        lo     = as.numeric(quantile(D2_reg, 0.1, na.rm = TRUE, names = FALSE)),
        hi     = as.numeric(quantile(D2_reg, 0.9, na.rm = TRUE, names = FALSE)),
        eps.lo = qs[1],
        eps.hi = qs[2]
    )

    ## where does eps_abc fall among the pairwise distances?
    eps_abc.pct <- if (!is.na(eps_abc)) mean(all.d < eps_abc) else NA_real_
    eps_abc.in.plateau <- !is.na(eps_abc) && eps_abc >= qs[1] && eps_abc <= qs[2]

    out <- list(
        eps_grid    = eps_grid,
        C           = C,
        C_reps      = C_reps,
        C_obs       = C_obs,
        eps_mid     = eps_mid,
        D2          = D2,
        D2_reps     = D2_reps,
        eps_abc     = eps_abc,
        eps_abc.pct = eps_abc.pct,
        eps_abc.in.plateau = eps_abc.in.plateau,
        plateau     = plateau,
        tol         = tol,
        distance    = distance,
        nsub        = m,
        nrep        = nrep,
        smooth      = smooth,
        n           = n.good,
        numstat     = numstat,
        names       = list(statistics.names = statnames),
        call        = match.call()
    )
    class(out) <- "corrdim4abc"
    invisible(out)
}


is.corrdim4abc <- function(x) inherits(x, "corrdim4abc")


plot.corrdim4abc <- function(x, file = NULL, postscript = FALSE, onefile = TRUE, ...) {

    if (!inherits(x, "corrdim4abc"))
        stop("Use only with objects of class \"corrdim4abc\".", call. = FALSE)

    ## Devices ('file' is a character basename, e.g. "corrdim").
    ## Open the target device BEFORE any other graphics call, otherwise the
    ## first plotting op leaks a blank Rplots.pdf to the default device in
    ## batch / Rscript runs.
    if (!is.null(file)) {
        file <- as.character(file)
        if (!postscript) pdf(file = paste(file, "pdf", sep = "."), onefile = onefile)
        else             postscript(file = paste(file, "ps", sep = "."), onefile = onefile)
        on.exit(dev.off(), add = TRUE)          # device (and its par) discarded on exit
        par(mfrow = c(1, 2), cex = 1, cex.main = 1.1, cex.lab = 1)
    } else {
        old_par <- par(mfrow = c(1, 2), cex = 1, cex.main = 1.1, cex.lab = 1)
        on.exit(par(old_par), add = TRUE)       # restore caller's device par
    }

    plat <- x$plateau

    ## -- Panel 1: local slope D2 vs eps (the headline figure) --------------
    ## clamp y to a sensible band: extrapolated slopes in the shaded tails
    ## can spike, so cap by a high quantile rather than the raw max.
    finiteD2 <- x$D2[is.finite(x$D2)]
    ymax <- max(plat$value * 1.5, plat$hi * 1.3,
                as.numeric(quantile(finiteD2, 0.95, names = FALSE)), na.rm = TRUE)
    ylim <- c(0, ymax)
    plot(x$eps_mid, x$D2, log = "x", type = "n", ylim = ylim,
         xlab = "epsilon  (ABC metric)", ylab = expression(D[2](epsilon)),
         main = "Local correlation dimension")

    ## shade the undersampled (small-eps) and saturated (large-eps) regimes
    usr <- par("usr")  # x in log10 units
    rect(10^usr[1], usr[3], plat$eps.lo, usr[4],
         col = adjustcolor("grey80", 0.4), border = NA)
    rect(plat$eps.hi, usr[3], 10^usr[2], usr[4],
         col = adjustcolor("grey80", 0.4), border = NA)

    ## faint per-replicate slopes
    if (x$nrep > 1)
        for (r in seq_len(x$nrep))
            lines(x$eps_mid, x$D2_reps[r, ], col = adjustcolor("steelblue", 0.25))

    lines(x$eps_mid, x$D2, col = "steelblue", lwd = 2)
    abline(h = plat$value, col = "black", lty = 2)

    if (!is.na(x$eps_abc)) {
        abline(v = x$eps_abc, col = "red", lwd = 2)
        mtext(expression(epsilon[abc]), side = 3, at = x$eps_abc, col = "red", line = 0.2)
    }
    legend("topright", bty = "n", cex = 0.8,
           legend = c(sprintf("plateau D2 = %.2f [%.2f, %.2f]",
                              plat$value, plat$lo, plat$hi),
                      "scaling region", "excluded (under/over)"),
           lty = c(2, NA, NA), pch = c(NA, 15, 15),
           col = c("black", "white", "grey80"))

    ## -- Panel 2: correlation integral log C vs log eps --------------------
    ## drop C == 0 points (small eps) so log axis doesn't warn / clip
    p <- x$C > 0
    plot(x$eps_grid[p], x$C[p], log = "xy", type = "n",
         xlab = "epsilon  (ABC metric)", ylab = "C(epsilon)",
         main = "Correlation integral")
    rect(plat$eps.lo, 10^par("usr")[3], plat$eps.hi, 10^par("usr")[4],
         col = adjustcolor("grey80", 0.3), border = NA)
    lines(x$eps_grid[p], x$C[p], col = "steelblue", lwd = 2)
    if (!is.null(x$C_obs)) {
        po <- x$C_obs > 0
        lines(x$eps_grid[po], x$C_obs[po], col = "darkgreen", lwd = 2, lty = 1)
    }
    if (!is.na(x$eps_abc))
        abline(v = x$eps_abc, col = "red", lwd = 2)
    legend("topleft", bty = "n", cex = 0.8,
           legend = c("all-pairs C(eps)",
                      if (!is.null(x$C_obs)) "pointwise C_obs(eps)  [ABC acceptance]"),
           col = c("steelblue", if (!is.null(x$C_obs)) "darkgreen"), lwd = 2)

    invisible()
}


summary.corrdim4abc <- function(object, print = TRUE,
                                digits = max(3, getOption("digits") - 3), ...) {

    if (!inherits(object, "corrdim4abc"))
        stop("Use only with objects of class \"corrdim4abc\".", call. = FALSE)

    x    <- object
    plat <- x$plateau

    if (print) {
        cat("Correlation dimension (D2) in ABC's metric\n")
        cat("-------------------------------------------\n")
        cat(sprintf("  distance            : %s\n", x$distance))
        cat(sprintf("  usable simulations  : %d   (summary stats: %d)\n", x$n, x$numstat))
        cat(sprintf("  subsample / reps    : %d points x %d replicate(s)\n", x$nsub, x$nrep))
        cat(sprintf("  scaling region (eps): [%.3g, %.3g]\n", plat$eps.lo, plat$eps.hi))
        cat("\n")
        cat(sprintf("  EFFECTIVE DIMENSION : %.2f   [%.2f, %.2f]\n",
                    round(plat$value, digits), round(plat$lo, digits), round(plat$hi, digits)))
        cat("    (noise-inclusive UPPER bound on the informative dimension:\n")
        cat(sprintf("     informative <= D2 <= nominal %d)\n", x$numstat))
        cat("\n")
        ## target is optional (off by default); only shown when supplied
        if (!is.na(x$eps_abc)) {
            cat("\n")
            cat(sprintf("  eps_abc (tol=%.3g)   : %.4g  (%.1f%% of pairwise distances)\n",
                        x$tol, x$eps_abc, 100 * x$eps_abc.pct))
            cat(sprintf("  eps_abc in plateau  : %s\n",
                        if (x$eps_abc.in.plateau) "YES  -- ABC operates at the D2 plateau"
                        else "NO   -- check sampling / scaling region"))
        }
    }

    res <- c(D2 = plat$value, lo = plat$lo, hi = plat$hi,
             eps_abc = x$eps_abc, eps_abc.pct = x$eps_abc.pct)
    invisible(res)
}
