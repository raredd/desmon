#' Compute Group Sequential Use Function Boundaries
#'
#' Computes the one-sided group sequential boundary for a specified use
#' function at specified analysis times.  Optionally, also computes an
#' asymmetric lower boundary based on repeated confidence interval (RCI)
#' monitoring for early stopping in favor of the null hypothesis.
#'
#' @details
#' Calculates the group sequential boundaries for repeated significance tests
#' in group sequential analysis of clinical trials.  The algorithm is based on
#' the use function approach proposed by Lan and DeMets (1983, Biometrika) and
#' investigated further by Kim and DeMets (1987, Biometrika).  The information
#' time corresponds to the proportion of statistical information, which is the
#' essentially the same as the number of failures in a proportional hazards
#' model for failure time endpoints.
#'
#' For the truncated O-F boundary (\code{use=6}), first the regular O-F
#' boundary is computed, but if the critical value is larger than the specified
#' truncation value, the truncated value is used instead.  The actual error
#' spent is computed, and the early over-spending is made up as quickly as
#' possible, after which the boundary is similar to the ordinary O-F boundary.
#'
#' With asymmetric monitoring, the study is stopped early in favor of the null
#' hypothesis if a `lower boundary' is crossed.  Here the lower boundary is
#' based on a repeated confidence interval on the log hazard ratio (see
#' Jennison and Turnbull, 1990).  The RCI is constructed using the critical
#' value from the one-sided boundary specified by \code{usel} and
#' \code{alphal}.  In the program, this is converted to a boundary on the
#' logrank statistic.  This requires information on the alternative and the
#' total planned information for the study, which is specified through the
#' parameter \code{eta}.  The value of \code{eta} can be obtained from the
#' functions \code{seqopr} and \code{lr.inf}.  Introducing a lower boundary
#' reduces the probability of crossing the upper boundary.  After determining
#' the lower boundary, the program computes the upper boundary taking into
#' account the lower boundary.  The upper boundary will thus be affected by the
#' parameters specified for the lower boundary.  It is quite easy to specify
#' incompatible combinations for the upper and lower boundary, especially if
#' \code{alpha} or \code{alphal} is very large or if \code{eta} is small.
#'
#' @param inf information times of analyses (length <= 30); must be positive,
#' increasing and <= 1
#' @param alpha one-sided significance level of the group sequential test
#' @param use the type of use function: 1=O'Brien-Fleming, 2=Pocock, 3=linear,
#' 4=one and a half, 5=quadratic, 6=truncated O'Brien-Fleming
#' @param eta The mean parameter on the Brownian motion process scale.  Only
#' needed for the RCI lower boundary.
#' @param alphal The one-sided significance level used in the RCI monitoring
#' for stopping in favor of the null.  The confidence level of the RCI is
#' \code{1-2*alphal}.  If \code{alphal <= 0}, only the upper boundary is
#' computed.
#' @param usel The use function for determining critical values for the RCI
#' lower boundary (same codes as \code{use})
#' @param oftr The significance level at which the truncated O-F boundary is
#' truncated (upper boundary)
#' @param oftrl The significance level at which the truncated O-F boundary is
#' truncated (RCI lower boundary)
#'
#' @return a matrix with \code{length(inf)} rows giving the critical values on
#' the standard normal scale at the specified information times.  If
#' \code{alphal<=0}, there is a single column giving the critical values for
#' the upper one-sided boundary.  If \code{alphal>0}, then there is a second
#' column giving the critical value for the lower boundary for early stopping
#' in favor of the null.  Note that these critical values are on the normalized
#' test statistic scale, which are NOT the critical values used in constructing
#' the RCI.
#' @note Interface to the Fortran code for the program \code{sequse}, which was
#' written by Kyungmann Kim and modified for truncated O-F boundaries and
#' asymmetric lower boundaries by Bob Gray
#'
#' @seealso
#' \code{\link{seqopr}}; \code{\link{lr.inf}}; \code{\link{seqp}}
#'
#' @references
#' Lan and DeMets (1983). \emph{Biometrika}.
#'
#' Kim and DeMets (1987). \emph{Biometrika}.
#'
#' Jennison and Turnbull (1990). \emph{Statistical Science} \strong{5}:299-317.
#'
#' @keywords design
#'
#' @examples
#' sequse((1:4) / 4)
#' sequse((1:4) / 4, use = 6)
#' sequse((1:4) / 4, use = 6, alphal = 0.025, eta = 2)
#'
#' @export sequse

sequse <- function(inf, alpha = 0.025, use = 6, eta = 0, alphal = 0,
                   usel = 6, oftr = alpha / 50, oftrl = alphal / 50) {
  # use=type use function 1=OF, 2=Pocock, 3=linear, 4=1^1.5, 5=t^2
  # 6=truncated O-F
  if (min(inf) <= 0 | max(inf) > 1)
    stop('inf must be >0, <=1')
  if (length(inf) > 30)
    stop('inf too long')
  if (alpha <= 0 | alpha >= 0.5)
    stop('alpha out of range')
  if (use < 1 | use > 6)
    stop('use must be >=1, <=6')
  inf <- sort(inf)
  z <- .Fortran(
    'squse', as.double(alpha), as.double(alphal), as.integer(length(inf)),
    as.double(c(0, inf)), as.integer(use), as.integer(usel), as.double(oftr),
    as.double(oftrl), as.double(eta), double(length(inf)), double(length(inf)),
    ierr = integer(1), PACKAGE = 'desmon'
  )[10:12]
  if (z$ierr > 0)
    stop('Boundaries could not be computed')
  if (alphal <= 0) {
    out <- as.matrix(z[[1L]])
    dimnames(out) <- list(format(round(inf, 2L)), 'upper')
    out
  } else {
    out <- cbind(z[[1L]], z[[2L]])
    dimnames(out) <- list(format(round(inf, 2L)), c('upper', 'lower'))
    out
  }
}
