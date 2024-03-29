% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intpow.R
\name{intpow}
\alias{intpow}
\alias{hazsf}
\alias{expinf}
\alias{siminf}
\title{Approximate Power of Cox Model Interaction Tests}
\usage{
intpow(
  n,
  prx,
  prz,
  ri,
  rt,
  rf,
  l0,
  acc.per,
  add.fu,
  alpha2 = 0.05,
  hazs = NULL
)
}
\arguments{
\item{n}{Sample size}

\item{prx}{Proportion randomized to treatment 1 (x=1)}

\item{prz}{Proportion positive for the marker (z=1)}

\item{ri}{Interaction ratio: the treatment hazard ratio (treatment 1 /
treatment 0) in the marker positive group divided by the treatment hazard
ratio in the marker negative group}

\item{rt}{Average treatment hazard ratio}

\item{rf}{Average marker hazard ratio (positive/negative)}

\item{l0}{Overall average hazard}

\item{acc.per}{Number of time units of accrual}

\item{add.fu}{Number of time units of follow-up after the end of accrual}

\item{alpha2}{Two-sided type I error rate for the interaction test}

\item{hazs}{Vector giving the constant hazard rates in the 4 groups, in the
order \code{(x,z)=(0,0),(0,1),(1,0),(1,1)}}
}
\value{
\code{intpow} returns a list with components: \item{power }{ The
power of the interaction test} \item{hazards }{ The vector of hazards from
\code{hazsf}} \item{var}{ The inverse expected information matrix (per
sample)} \item{nevents}{ The expected number of events in each of the
\code{x} by \code{z} combinations}

\code{hazsf} returns the vector of hazards meeting the constraints on the
ratios and average hazard, in the order \code{(x,z)=(0,0),(0,1),(1,0),(1,1)}

\code{expinf} returns the expected information matrix per subject

\code{siminf} returns a vector containing the Wald test statistic and the
estimated variance of the interaction term
}
\description{
Approximates the power of the large sample partial likelihood tests for a
treatment by biologic marker interaction test in a Cox model using the
normal approximation, assuming exponential distributions.
}
\details{
For a Cox model with two binary factors (x=treatment and z=biologic marker),
\code{intpow} computes the power of the large sample partial likelihood
tests for no interaction. Constant hazards are assumed within each of the 4
covariate combinations.

Let \code{p(x,z)} be the sample proportion and \code{h(x,z)} be the hazard
for treatment \code{x} in marker group \code{z}. Treatment is assumed
independent of the marker, so the proportion of the sample in the 4
combinations is \code{p(0,0)=(1-prx)*(1-prz),
p(0,1)=(1-prx)*prz, p(1,0)=prx*(1-prz), p(1,1)=prx*prz}.

The hazards can either be given explicitly (argument \code{hazs}) or
through the ratios and the overall average hazard. The ratios are defined by
\code{ri=h(1,1)*h(0,0)/(h(0,1)*h(1,0))},
\code{rt=((1-prz)*h(1,0)+prz*h(1,1))/((1-prz)*h(0,0)+prz*h(0,1))}, and
\code{rf=((1-prx)*h(0,1)+prx*h(1,1))/((1-prx)*h(0,0)+prx*h(1,0))}.

Given a set of ratios and an average hazard \code{l0}, \code{hazsf} computes
the hazards in the individual groups to give the specified ratios and to
satisfy the constraint
\code{p(0,0)*h(0,0)+p(0,1)*h(0,1)+p(1,0)*h(1,0)+p(1,1)*h(1,1)=l0}.
\code{hazsf} is called by \code{intpow} if \code{hazs} is not specified.

The calculations assume a clinical trial setting with accrual uniform over
the period \code{(0,acc.per)} and an additional \code{add.fu} units of
follow-up after completion of accrual. Thus censoring is assumed uniform on
\code{(add.fu, add.fu+acc.per)}.

\code{expinf} computes the expected information matrix (per subject) from
the Cox partial likelihood containing the covariates \code{x}, \code{z},
and \code{x*z}.

\code{siminf} calculates the Wald statistic for the interaction test and
the partial likelihood variance estimate of the variance of the interaction
coefficient for a single simulated sample from this model.

\code{intpow} calculates the power for the Wald test for the null hypothesis
that coefficient of the interaction term = 0 using the large sample normal
approximation. In large samples, this test is equivalent to the score and
partial likelihood ratio tests for this hypothesis.
}
\examples{
intpow(n = 2800, prx = 0.5, prz = 0.25, ri = 0.6, rt = 0.8, rf = 2,
       l0 = 0.045, acc.per = 2.5, add.fu = 3)

}
\keyword{design}
\keyword{survival}
