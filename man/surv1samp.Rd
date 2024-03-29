% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv1samp.R
\name{surv1samp}
\alias{surv1samp}
\title{Sample Size Calculation for a One-Sample Exponential}
\usage{
surv1samp(
  control.rate = NULL,
  pct.imp = NULL,
  acc.per,
  add.fu,
  alpha = 0.1,
  power = 0.9,
  control.med = NULL,
  control.survp = NULL,
  tp = NULL,
  pct.reduc = NULL,
  nonpar = FALSE,
  alt.rate = NULL,
  alt.med = NULL,
  alt.survp = NULL
)
}
\arguments{
\item{control.rate}{The hazard rate under the null hypothesis}

\item{pct.imp}{The percent improvement in the median failure time under the
alternative}

\item{acc.per}{The length of the time period during which patients are
accrued}

\item{add.fu}{The length of time after the end of accrual during which
patients continue to be followed prior to analysis}

\item{alpha}{The one-sided type I error rate}

\item{power}{The power of the test}

\item{control.med}{The median under the null}

\item{control.survp}{The survival probability at \code{tp} under the null}

\item{tp}{Time point for survival probability values and/or the
nonparametric test}

\item{pct.reduc}{Percent reduction in the hazard rate under the alternative}

\item{nonpar}{If TRUE, the test based on the nonparametric cumulative hazard
estimator at \code{tp} is used}

\item{alt.rate}{The hazard rate under the alternative}

\item{alt.med}{The median under the alternative}

\item{alt.survp}{The survival probability at \code{tp} under the alternative}
}
\value{
A vector giving the sample size (n), the number of failures needed
(nd), the proportion of failures (fail.prob), the hazard ratio of the null
to alternative hazard rates (haz.ratio), the input values of \code{pct.imp}
and \code{control.rate}, the hazard rate under the alternative (alt.rate)
and the input values of \code{alpha} and \code{beta}.
}
\description{
Computes the sample size for a given size and power for the test for a
specified improvement in the exponential rate parameter.
}
\details{
Accrual is assumed to be uniform over the accrual period, so censoring is
uniform on the interval \code{add.fu} to \code{acc.per + add.fu}. The power
is based on the Wald test for the log failure rate parameter unless
\code{nonpar=TRUE}, in which case the test is based on the nonparametric
estimate of the cumulative hazard at \code{tp}. Simulations suggest that the
distribution approximation to the latter test is anti-conservative, so the
test (and the corresponding power calculations should be used with extreme
caution). Exponential failures are assumed for the calculations, even for
the nonparametric test.

Any combination of parameters determining the null and alternative event
rates can be specified.
}
\examples{
surv1samp(-1.333 * log(0.96) / 5, 33.3, 4, 3, alpha = 0.1, power = 0.95)
surv1samp(-log(0.95) / 5, 100 * (log(0.95) / log(0.96) - 1), 4, 3,
          alpha = 0.1, power = 0.95)
surv1samp(-log(0.955) / 5, 100 * (log(0.955) / log(0.965) - 1), 4, 3,
          alpha = 0.1, power = 0.95)

}
\keyword{design}
\keyword{survival}
