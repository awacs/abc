\name{corrDim4abc}
\alias{corrDim4abc}
\alias{is.corrdim4abc}
\alias{plot.corrdim4abc}
\alias{summary.corrdim4abc}
\title{
  Correlation dimension (D2) of the simulated summary statistics in ABC's metric
}
\description{
  Estimates the Grassberger-Procaccia correlation dimension (\eqn{D_2}) of the
  simulated summary-statistic cloud, computed in the exact representation and
  distance metric used by \code{\link{abc}}. ABC's acceptance count within a
  tolerance \eqn{\epsilon} is the correlation integral \eqn{C(\epsilon)}, whose
  scaling exponent is \eqn{D_2}; the reported value is therefore the
  effective dimension ABC operates in at its acceptance scale, and a
  noise-inclusive \emph{upper} bound on the data's informative dimension
  (informative \eqn{\le D_2 \le} nominal number of statistics).
}
\usage{
corrDim4abc(sumstat, target = NULL, distance = "euclidean", subset = NULL,
            tol = 0.01, nsub = 4000, nrep = 1, ngrid = 250, smooth = TRUE,
            scaling.quantiles = c(0.05, 0.40), seed = NULL, \dots)
}
\arguments{
  \item{sumstat}{a matrix or data frame of simulated summary statistics, in the
    same representation passed to \code{\link{abc}} (e.g. already normalised/cut).}
  \item{target}{optional observed summary statistics, \code{NULL} (off) by
    default. The correlation dimension is a property of the simulated cloud
    alone; supplying \code{target} only adds the operating-point marker
    \code{eps_abc} and the pointwise acceptance curve \eqn{C_{obs}(\epsilon)}.}
  \item{distance}{distance metric, matching the \code{abc()} run: one of
    \code{"euclidean"}, \code{"manhattan"} or \code{"wasserstein"}. The same
    per-column MAD weighting / CDF transform as \code{\link{abc}} is applied.}
  \item{subset}{optional logical vector to restrict the simulations used.}
  \item{tol}{ABC acceptance rate used to locate \code{eps_abc} (default 0.01).}
  \item{nsub}{number of points subsampled for the all-pairs correlation
    integral (cost is quadratic in \code{nsub}).}
  \item{nrep}{number of independent subsamples, for a robustness band.}
  \item{ngrid}{number of log-spaced \eqn{\epsilon} values.}
  \item{smooth}{if \code{TRUE}, the local slope is the derivative of a smoothing
    spline of \eqn{\log C} vs \eqn{\log\epsilon}; otherwise finite differences.}
  \item{scaling.quantiles}{lower/upper quantiles of the pairwise distances that
    delimit the plateau (scaling) region; smaller/larger scales are excluded.}
  \item{seed}{optional random seed.}
  \item{\dots}{unused.}
}
\value{
  An object of class \code{"corrdim4abc"}: a list with the \eqn{\epsilon} grid,
  the correlation integral \code{C} (pooled and per replicate), the optional
  pointwise curve \code{C_obs}, the local slope \code{D2}, the operating point
  \code{eps_abc}, and a \code{plateau} summary (median height with a range).
  Methods \code{plot} and \code{summary} are provided.
}
\details{
  The metric is reproduced from the package internals (\code{calc_distance},
  \code{normalise}): \code{"euclidean"}/\code{"manhattan"} apply per-column
  bank-MAD weighting, while \code{"wasserstein"} uses per-row CDFs without
  weighting. \code{eps_abc} follows the same tolerance-to-quantile rule as
  \code{\link{abc}}.
}
\seealso{
  \code{\link{abc}}, \code{\link{cv4abc}}
}
\examples{
\dontrun{
## dimension of the simulated cloud (no target needed)
cd <- corrDim4abc(sumstat = normalized_csfs_sim,
                  distance = "euclidean", nsub = 4000, nrep = 3)
summary(cd)
plot(cd, file = "corrdim")

## optional: also mark the ABC operating point eps_abc
cd2 <- corrDim4abc(sumstat = normalized_csfs_sim, target = tgt,
                   distance = "euclidean", tol = 0.01, nsub = 4000, nrep = 3)
}
}
\keyword{models}
