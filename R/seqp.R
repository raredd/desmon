#' Compute Boundary Crossing Probabilities
#'
#' Computes probabilities of crossing arbitrary user specified upper and/or
#' lower boundaries
#'
#' @details
#' For a group sequential design with arbitrary boundaries, computes the
#' probability of crossing the upper and lower boundary at each analysis.  The
#' probability at each analysis is the probability of not crossing either
#' boundary at earlier analyses and crossing at the current analysis.
#'
#' \code{eta} should be 0 to compute probabilities under the null, and should be
#' set to a value which is a function of the alternative on a transformed scale
#' to compute probabilities under the alternative.  For survival studies, eta
#' can be computed using the functions \code{seqopr} or \code{lr.inf}.
#' (\code{eta} is the mean of Y(1), where Y(t) is the partial sum process
#' standardized to have variance 1 at t=1.)
#'
#' @param inf.times Information times of analyses
#' @param upper Upper boundary (for rejecting H0) on standard normal scale
#' @param lower Lower boundary (for rejecting HA) on standard normal scale
#' @param eta The mean parameter on the Brownian motion process scale.
#'
#' @return
#' A \code{length(inf.times)} by 5 matrix with columns giving the
#' information times, the upper boundary, the lower boundary, the probabilities
#' of crossing the upper boundary, and the probabilities of crossing the lower
#' boundary.
#'
#' @seealso
#' \code{\link{sequse}}; \code{\link{seqopr}}; \code{\link{lr.inf}}
#'
#' @keywords design
#'
#' @examples
#' seqp((1:4) / 4, c(2, 2, 3, 2), eta = 3)
#' seqp((1:4) / 4, c(2, 2, 3, 2), eta = 0)
#'
#' @export seqp

seqp <- function(inf.times, upper, lower = NULL, eta = 0) {
  # compute boundary crossing prob for arbitrary sequential boundary
  kk <- length(inf.times)
  upper <- if (is.null(upper))
    rep(6, kk) else pmin(6, upper)
  lower <- if (is.null(lower))
    rep(-6, kk) else pmax(lower, -6)
  if (kk > 30)
    stop('too many analyses')
  if (length(upper) != kk | length(lower) != kk)
    stop('boundaries wrong length')

  u2 <- .Fortran(
    'sqopr3', kk = as.integer(kk), pt = as.double(c(0, inf.times)),
    eta = as.double(eta), uz = as.double(upper), uzl = as.double(lower),
    pnu = double(kk + 1), pnl = double(kk + 1), ierr = integer(1),
    PACKAGE = 'desmon'
  )
  if (u2$ierr > 0)
    stop('Boundaries could not be computed')
  w2 <- cbind(u2$pt[-1L], u2$uz, u2$uzl, diff(u2$pnu), diff(u2$pnl))
  dimnames(w2) <- list(NULL, c('inf.times', 'Uz', 'Lz', 'P.ge.Uz', 'P.le.Lz'))
  w2
}
