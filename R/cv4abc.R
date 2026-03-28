######################################################################
#
# cv4abc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum
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
# Contains: cv4abc, is.cv4abc, summary.cv4abc, plot.cv4abc
#
######################################################################


## SBC coverage levels, matching validate_npe.py
.SBC_LEVELS <- round(seq(0.05, 0.95, length.out = 19), 3)

## Internal helper: coverage matrix from rank matrix
## ranks : (nval x np) matrix of normalised ranks in [0,1]
## returns (nlevels x np) matrix of actual coverage at each .SBC_LEVELS level
.sbc_coverage <- function(ranks, levels = .SBC_LEVELS) {
    np       <- ncol(ranks)
    coverage <- matrix(0, nrow = length(levels), ncol = np)
    for (li in seq_along(levels)) {
        lo <- (1 - levels[li]) / 2
        hi <- 1 - lo
        coverage[li, ] <- colMeans(ranks >= lo & ranks <= hi, na.rm = TRUE)
    }
    coverage
}


cv4abc <- function(param, sumstat, abc.out = NULL, nval, tols,
                   statistic = "median", prior.range = NULL,
                   method, hcorr = TRUE,
                   transf = "none", logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov",
                   distance = "euclidean", numnet = 10, sizenet = 5,
                   lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 500,
                   keep_posterior = FALSE,
                   ncores = 1,
                   ...){

    mywarn <- options()$warn
    options(warn=-1)
    linout <- TRUE
    ## checks:
    ## ########
    if(!any(statistic == c("median", "mean", "mode"))){
        stop("Statistic has to be mean, median or mode.", call.=F)
    }
    if(is.null(abc.out) && missing(method)) stop("Method must be supplied when 'abc.out' is NULL.", call.=F)
    if(missing(nval)) stop("'nval' must be supplied.", call.=F)

    ## set random seeds
    ## ################
    if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)

    ## define defaults:
    ## #################

    if(!is.null(abc.out)){
        subset       <- abc.out$na.action
        method       <- abc.out$method
        transf       <- abc.out$transf
        logit.bounds <- abc.out$logit.bounds
        kernel       <- "epanechnikov"
        distance     <- abc.out$distance
    }

    ## checks, numbers of stats, params, sims
    if(is.null(dim(param))){
        np    <- 1
        param <- as.data.frame(param)
        names(param) <- "P1"
    }
    else np <- dim(param)[2]
    if(is.null(dim(sumstat))){
        numstat <- 1
        sumstat <- as.data.frame(sumstat)
        names(sumstat) <- "S1"
    }
    else{
        numstat <- dim(sumstat)[2]
    }
    numsim <- dim(sumstat)[1]

    ## paramnames & statnames
    if(!is.null(abc.out)){ # paramnames & statnames from abc.out
        if(np != abc.out$numparam || numstat != abc.out$numstat || numsim != length(abc.out$na.action)){
            stop("The number of parameters, summary statistics, or simulations provided in 'param' or 'sumstat' are not the same as in 'abc.out'.", call.=F)
        }
        else if(!prod(colnames(param) %in% abc.out$names$parameter.names)){
            stop("Parameters in 'param' are not the same as in 'abc.out', or different names are used.", call.=F)
        }
        else if(!prod(colnames(sumstat) %in% abc.out$names$statistics.names)){
            stop("Summary statistics in 'sumstat' are not the same as in 'abc.out', or different names are used.", call.=F)
        }
        else{
            paramnames <- abc.out$names$parameter.names
            statnames  <- abc.out$names$statistics.names
        }
    }
    else{ # give paramnames & statnames o/w
        if(length(colnames(param))){
            paramnames <- colnames(param)
        }
        else{
            paramnames <- paste("P", 1:np, sep="")
        }
    }
    if(length(colnames(sumstat))){
        statnames <- colnames(sumstat)
    }
    else{
        statnames <- paste("S", 1:numstat, sep="")
    }

    ## indexes for the CV sample and check that the sample is not actually an NA
    gwt <- rep(TRUE, length(sumstat[,1]))
    gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
    if(is.null(subset)) subset <- rep(TRUE, length(sumstat[,1]))
    gwt    <- as.logical(gwt * subset)
    cvsamp <- sample(1:numsim, nval, prob = gwt/sum(gwt))

    tols      <- sort(tols)
    num.panel <- length(tols)
    alltol    <- list()   # point estimates  (nval x np), keyed by tolerance
    allranks  <- list()   # SBC ranks        (nval x np), keyed by tolerance
    mycall    <- list()
    allpost   <- list()   # raw posteriors, only populated when keep_posterior = TRUE

    for(mytol in tols){
        res_estim <- matrix(ncol = np, nrow = nval)
        res_ranks <- matrix(ncol = np, nrow = nval)
        post_list <- if (keep_posterior) vector("list", nval) else NULL

        cv_one <- function(i) {
            mysamp    <- cvsamp[i]
            mytrue    <- param[mysamp, ]
            mytarget  <- sumstat[mysamp, ]
            myparam   <- param[-mysamp, ]
            mysumstat <- sumstat[-mysamp, ]
            mysubset  <- subset[-mysamp]

            subres <- withCallingHandlers(
                abc(target = mytarget, param = myparam, sumstat = mysumstat, tol = mytol,
                    subset = mysubset, method = method, transf = transf,
                    logit.bounds = logit.bounds, kernel = kernel, hcorr = hcorr,
                    distance = distance),
                warning = namesWarningFilter)

            ## --- point estimate (original behaviour) ---
            if(statistic == "median") estim <- invisible(summary.abc(subres, print = F, ...)[3, ])
            if(statistic == "mean")   estim <- invisible(summary.abc(subres, print = F, ...)[4, ])
            if(statistic == "mode")   estim <- invisible(summary.abc(subres, print = F, ...)[5, ])

            ## --- SBC rank ---
            post <- if (!is.null(subres$adj.values)) subres$adj.values else subres$unadj.values
            post <- matrix(post, ncol = np)
            true_row <- as.numeric(mytrue)
            rank_i <- colMeans(post < matrix(true_row, nrow = nrow(post), ncol = np, byrow = TRUE))

            list(estim = estim, rank = rank_i,
                 post = if (keep_posterior) post else NULL,
                 subres = subres)
        }

        if (ncores > 1) {
            results <- parallel::mclapply(1:nval, cv_one, mc.cores = ncores)
        } else {
            results <- lapply(1:nval, cv_one)
        }

        for (i in 1:nval) {
            res_estim[i, ] <- results[[i]]$estim
            res_ranks[i, ] <- results[[i]]$rank
            if (keep_posterior) post_list[[i]] <- results[[i]]$post
        }
        subres <- results[[nval]]$subres

        if(np == 1){
            res_estim <- c(res_estim)
            res_ranks <- c(res_ranks)
        } else {
            colnames(res_estim) <- paramnames
            colnames(res_ranks) <- paramnames
        }

        tkey <- paste("tol", mytol, sep="")
        alltol[[tkey]]   <- res_estim
        allranks[[tkey]] <- res_ranks
        mycall[[tkey]]   <- call("abc",
                                 target = quote(target), param = quote(param), sumstat = quote(sumstat),
                                 tol = mytol, subset = quote(subset),
                                 method = subres$method, hcorr = subres$hcorr, transf = subres$transf,
                                 logit.bounds = subres$logit.bounds, kernel = subres$kernel)
        if (keep_posterior) allpost[[tkey]] <- post_list
    }

    if(np == 1){
        true <- as.data.frame(param[cvsamp, ])
        names(true) <- paramnames
    }
    else true <- param[cvsamp, ]

    cv4abc.out <- list(
        calls      = mycall,
        cvsamples  = cvsamp,
        tols       = tols,
        true       = true,
        estim      = alltol,
        ranks      = allranks,
        posteriors = if (keep_posterior) allpost else NULL,
        names      = list(parameter.names = paramnames, statistics.names = statnames),
        seed       = seed
    )

    options(warn = mywarn)
    class(cv4abc.out) <- "cv4abc"
    invisible(cv4abc.out)
}


is.cv4abc <- function(x){
    if (inherits(x, "cv4abc")) TRUE
    else FALSE
}


plot.cv4abc <- function(x, sbc = FALSE, exclude = NULL, log = NULL, file = NULL,
                        postscript = FALSE, onefile = TRUE,
                        ask = !is.null(deviceIsInteractive()), caption = NULL, bins = 20, ...){

    if (!inherits(x, "cv4abc"))
        stop("Use only with objects of class \"cv4abc\".", call.=F)

    cv4abc.out <- x
    tols     <- cv4abc.out$tols
    numtols  <- length(tols)
    np       <- length(cv4abc.out$names$parameter.names)
    parnames <- cv4abc.out$names$parameter.names
    nval     <- length(cv4abc.out$cvsamples)

    if(is.null(caption)) caption <- as.graphicsAnnot(parnames)

    ## Devices
    save.devAskNewPage <- devAskNewPage()
    if(!is.null(file)){
        file <- substitute(file)
        if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile = onefile)
        if(postscript)  postscript(file = paste(file, "ps",  sep="."), onefile = onefile)
    }

    if(sbc){
        ## -- SBC mode: rank histograms + coverage curve -------------------------
        levels <- .SBC_LEVELS
        cols   <- rainbow(np)

        for(mytol in tols){
            tkey  <- paste("tol", mytol, sep="")
            ranks <- matrix(cv4abc.out$ranks[[tkey]], ncol = np)

            if (!is.null(file) || !ask) {
                ## rank histograms
                old_par <- par(mfrow = c(np, 1), cex = 1, cex.main = 1.1, cex.lab = 1)
                on.exit(par(old_par), add = TRUE)
            } else {
                devAskNewPage(TRUE)
                old_par <- par(mfrow = c(np, 1), cex = 1, cex.main = 1.1, cex.lab = 1)
                on.exit(par(old_par), add = TRUE)
            }

            expected_height <- nval / bins
            for(j in 1:np){
                hist(ranks[, j], breaks = seq(0, 1, length.out = bins + 1),
                     col = "steelblue", border = "white",
                     main = paste0(caption[j], "\n(tol=", mytol, ")"),
                     xlab = "normalised rank", ylab = if(j==1) "count" else "")
                abline(h = expected_height, col = "red", lty = 2, lwd = 1.5)
            }
            mtext("SBC Rank Histograms  (flat = well calibrated)",
                  outer = TRUE, line = -1.5, cex = 1.1)

            ## coverage curve
            coverage <- .sbc_coverage(ranks, levels)
            par(mfrow = c(1, 1))
            plot(NULL, xlim = c(0,1), ylim = c(0,1),
                 xlab = "Expected coverage", ylab = "Actual coverage",
                 main = paste0("SBC Coverage Curve  (tol=", mytol, ")\n",
                               "closer to diagonal = better calibrated"))
            abline(0, 1, lty = 2, col = "black", lwd = 1)
            for(j in 1:np)
                lines(levels, coverage[, j], col = cols[j], lwd = 2, type = "b", pch = 16, cex = 0.7)
            legend("topleft", legend = parnames, col = cols, lwd = 2, cex = 0.8)
        }

    } else {
        ## -- Normal mode: true vs estimated scatter ------------------------------
        true  <- cv4abc.out$true
        estim <- as.data.frame(cv4abc.out$estim)

        if(is.null(log)) log <- rep("", each = np)
        else if(length(log) != np) stop("error in argument 'log': provide scale for all parameters.", call.=F)

        cols <- heat.colors(numtols)
        pch  <- 20

        if(is.null(file) && ask && 1 < np) devAskNewPage(TRUE)

        par(cex = 1, cex.main = 1.2, cex.lab = 1.1)
        for(par in 1:np){
            mylog   <- log[par]
            columns <- seq(par, numtols * np, by = np)
            if(!is.null(exclude)){
                plot(rep(true[-exclude, par], times = numtols),
                     unlist(estim[-exclude, columns]),
                     col = rep(cols, each = nval - length(exclude)),
                     pch = pch, log = mylog,
                     xlab = "True value", ylab = "Estimated value", main = caption[par])
            } else {
                plot(rep(true[, par], times = numtols),
                     unlist(estim[, columns]),
                     col = rep(cols, each = nval),
                     pch = pch, log = mylog,
                     xlab = "True value", ylab = "Estimated value", main = caption[par])
            }
            abline(0, 1)
        }
    }

    if(!is.null(file)) dev.off()
    else devAskNewPage(save.devAskNewPage)

    invisible()
}


summary.cv4abc <- function(object, sbc = FALSE, print = TRUE,
                           digits = max(3, getOption("digits")-3), ...){

    if (!inherits(object, "cv4abc"))
        stop("Use only with objects of class \"cv4abc\".", call.=F)

    cv4abc.out <- object
    tols     <- cv4abc.out$tols
    np       <- length(cv4abc.out$names$parameter.names)
    parnames <- cv4abc.out$names$parameter.names
    nval     <- length(cv4abc.out$cvsamples)

    if(sbc){
        ## -- SBC mode: coverage table + MAE from diagonal -----------------------
        levels <- .SBC_LEVELS
        cat(paste0("SBC calibration check based on ", nval, " cross-validation samples\n"))

        for(mytol in tols){
            tkey     <- paste("tol", mytol, sep="")
            ranks    <- matrix(cv4abc.out$ranks[[tkey]], ncol = np)
            coverage <- .sbc_coverage(ranks, levels)

            cat(paste0("\n-- Tolerance: ", mytol, " --\n"))
            sep <- paste(rep("-", 10 + 10 * np), collapse="")
            cat(sep, "\n")
            cat(sprintf("  %-8s", "Level"), sprintf("%10s", parnames), "\n")
            cat(sep, "\n")
            for(li in seq_along(levels)){
                cat(sprintf("  %-8s", paste0(round(levels[li]*100), "%")),
                    sprintf("%10.3f", coverage[li, ]), "\n")
            }
            cat(sep, "\n")

            expected <- matrix(levels, nrow = length(levels), ncol = np)
            mae      <- colMeans(abs(coverage - expected), na.rm = TRUE)
            cat("\n  Mean |actual - expected| coverage (lower = better calibrated):\n")
            for(j in seq_len(np)){
                flag <- if(mae[j] < 0.05) "well calibrated" else if(mae[j] < 0.10) "moderate" else "POOR"
                cat(sprintf("    %-14s %.4f  (%s)\n", parnames[j], mae[j], flag))
            }
        }
        invisible(cv4abc.out)

    } else {
        ## -- Normal mode: prediction error table (original behaviour) ------------
        true  <- cv4abc.out$true
        estim <- cv4abc.out$estim

        cat(paste("Prediction error based on a cross-validation sample of ", nval, "\n\n", sep=""))

        sqdiff  <- lapply(estim, function(a) apply((a - true)^2, 2, sum))
        truevar <- apply(true, 2, var) * nval
        prederr <- lapply(sqdiff, function(a) a / truevar)
        prederr <- t(as.data.frame(prederr))
        rownames(prederr) <- tols

        class(prederr) <- "table"
        prederr
    }
}
