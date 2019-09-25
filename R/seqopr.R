#' Compute Operating Characteristics for Group Sequential Designs
#' 
#' Computes operating characteristics (seqopr) and sample size (seqss) for
#' group sequential logrank tests in the setting of phase III studies with
#' failure time endpoints.
#' 
#' @aliases seqopr seqss lr.inf
#' 
#' @details
#' Assumes accrual is uniform at the rate \code{acc.rate} for \code{acc.per}
#' units of time, so the total number of subjects is \code{acc.rate * acc.per},
#' with interim analyses to be performed at the times in \code{a.times}.  The
#' final analysis is assumed to take place at the final time in \code{a.times}.
#' Significance levels at the analyses are determined from the error spending
#' rate function approach (see \code{sequse} for more details).  Optionally, an
#' asymmetric lower boundary for early stopping in favor of the null can be
#' included by specifying a value of \code{alphal > 0}.  An approximate lower
#' boundary based on repeated confidence intervals is used (see \code{sequse}).
#' Failure times are assumed to be exponentially distributed, with the
#' exponential failure hazard rates specified through the \code{control.rate}
#' parameter, or alternately the \code{control.med} parameter.  A stratified
#' sample with different failure rates in the strata can be specified by giving
#' a vector of control group hazard rates corresponding to the different
#' strata, and specifying the stratum proportions in \code{prop.strata}.  These
#' functions assume proportional hazards alternatives.  The improvement from
#' using the experimental treatment instead of control is specified as the
#' PERCENT improvement in the median time to failure.  This is related to the
#' hazard ratio through hazard ratio = \code{1/(1+pct.imp/100)}.  Alternately,
#' the percent reduction in the hazard rate (\code{pct.reduc}) can be
#' specified.
#' 
#' Given the specifications above, \code{seqopr} computes the boundaries,
#' power, boundary crossing probabilities and expected stopping times.  The
#' expected stopping times under the null are computed using the specified
#' calendar times of analysis and the null failure rates, but using the
#' information times and boundaries based on the expected information at these
#' times as determined under the alternative.  They therefore do not typically
#' give the correct expected stopping times under the design as it will be
#' implemented in practice.
#' 
#' In \code{seqss}, the desired power is specified, and a range is given for
#' the accrual period.  The program then searches for the accrual period that
#' gives the desired power.  The timing of the final analysis after the
#' completion of accrual (parameter \code{add.fu}) is kept fixed.  Since the
#' study duration is variable, an approximate interim analysis schedule is
#' specified through the parameters \code{anal.int} and \code{first.anal}.
#' When early stopping in favor of the null is included, it is possible that
#' some designs examined in the search will have incompatible specifications.
#' If this happens, try reducing the upper limit on \code{acc.per}.  If the
#' range of \code{acc.per} specified does not include a solution, then the
#' algorithm should give an error message stating that the values at the end
#' points are not of opposite signs.  If this occurs, try increasing the width
#' of the interval.
#' 
#' Given accrual and follow-up information and failure rates, \code{lr.inf}
#' computes the total planned information and the Brownian motion mean
#' parameter \code{eta}.  The main use of this function is to get the value of
#' \code{eta} for computing the lower boundary in \code{sequse}.
#' 
#' @usage
#' seqopr(acc.per, acc.rate, control.rate = NULL, pct.imp = NULL, 
#'        a.times, prop.strat = rep(1, length(control.rate)), alpha = 0.025,
#'        alphal = 0, p.con = 0.5, use = 6, usel = 6, oftr = alpha/50,
#'        oftrl = alphal/50, control.med = NULL, pct.reduc = NULL)
#'        
#' seqss(power, acc.per, acc.rate, add.fu, anal.int, first.anal = anal.int,
#'       control.rate = NULL, pct.imp = NULL, prop.strat = rep(1, length(control.rate)),
#'       alpha = 0.025, alphal = 0, p.con = 0.5, use = 6, usel = 6,
#'       oftr = alpha/50, oftrl = alphal/50, control.med = NULL, pct.reduc = NULL)
#' 
#' lr.inf(acc.per, acc.rate, add.fu, control.rate = NULL, pct.imp = NULL,
#'        prop.strat = rep(1, length(control.rate)), p.con = 0.5,
#'        control.med = NULL, pct.reduc = NULL)
#' 
#' @param acc.per In \code{sequse} and \code{lr.inf}, the planned accrual
#' period for the study.  In \code{seqss}, a vector of length 2 giving the min
#' and max values of the range of possible accrual periods over which to
#' search.
#' @param acc.rate Expected accrual rate
#' @param control.rate Vector giving the hazard rate on the control treatment
#' in each stratum
#' @param pct.imp The percent improvement in the median failure time for the
#' experimental treatment over the control
#' @param a.times Planned chronological times of interim analyses.  Must be
#' positive and in increasing order.
#' @param prop.strat Vector giving the proportion of the total sample in each
#' stratum.  The vector is renormalized to sum to 1, so the values only need to
#' be proportional to the proportions.
#' @param alpha The one-sided significance level of the group sequential test
#' @param alphal The one-sided significance level used in the RCI monitoring
#' for stopping in favor of the null.  The confidence level of the RCI is
#' \code{1-2*alphal}.  If \code{alphal <= 0}, only the upper boundary is used.
#' @param p.con The proportion of subjects to be randomized to the control
#' group
#' @param use the type of use function: 1=O'Brien-Fleming, 2=Pocock, 3=linear,
#' 4=one and a half, 5=quadratic, 6=truncated O'Brien-Fleming
#' @param usel The use function for determining critical values for the RCI
#' lower boundary (same codes as \code{use})
#' @param oftr The significance level at which the truncated O-F boundary is
#' truncated (upper boundary)
#' @param oftrl The significance level at which the truncated O-F boundary is
#' truncated (RCI lower boundary)
#' @param control.med Vector giving the median failure times on the control
#' treatment in each stratum (may be given instead of \code{control.rate})
#' @param pct.reduc The percent reduction in the hazard rate for the
#' experimental treatment relative to the control (may be given instead of
#' \code{pct.imp})
#' @param power The desired power for the study
#' @param add.fu The additional follow-up planned after the end of accrual
#' before the final analysis will be conducted
#' @param anal.int Interim analyses are planned to take place every
#' \code{anal.int} time units on the chronological time scale, beginning at
#' \code{first.anal}
#' @param first.anal Time of the first planned interim analysis
#' 
#' @return 
#' \code{seqopr} returns a list with the following components (except
#' \code{acc.per}).  \code{seqss} returns a list with the following components
#' computed from the final design.  \item{Expected.inf }{A matrix giving the
#' expected information times at the planned analysis times and the expected
#' number of failures under the null and alternative and the sample size at
#' each analysis.} \item{Boundaries }{A matrix giving the boundaries on the
#' standard normal test statistic scale at the expected information at the
#' planned analysis times, plus columns giving the upper and lower boundary
#' crossing probabilities under H0 and H1 at each analysis} \item{Rej.Probs }{A
#' vector of 4 values giving the size and power (upper boundary crossing
#' probabilities under H0 and H1) and the total lower boundary crossing
#' probabilities under H0 and H1. } \item{stop.times }{A 2 x 4 matrix of
#' expected stopping times} \item{eta }{The Brownian motion mean parameter}
#' \item{call }{The call to \code{seqopr}} \item{acc.per}{The accrual period
#' that gives the desired power}
#' 
#' \code{lr.inf} returns a vector giving the values of the Brownian motion mean
#' parameter \code{eta} and the expected number of failures under the
#' alternative.
#' 
#' @seealso \code{\link{sequse}}, \code{\link{powlgrnk6}}, \code{\link{seqp}}
#' @keywords design
#' 
#' @examples
#' seqopr(3,200,.1,50,c(2,3,4,5))
#' seqss(.8,c(2,6),200,2,.5,2,.1,50)
#' lr.inf(3,200,2,.1,50)
#' 
#' @export seqopr

seqopr <- function(acc.per,acc.rate,control.rate=NULL,pct.imp=NULL,a.times,
                   prop.strat=rep(1,length(control.rate)),
                   alpha=.025,alphal=0,p.con=.5,use=6,
                   usel=6,oftr=alpha/50,oftrl=alphal/50,
                   control.med=NULL,pct.reduc=NULL) {
  if (is.null(control.rate)) {
    if (is.null(control.med)) stop('control.rate or control.med must be given')
    control.rate <- log(2)/control.med
  }
  if (is.null(pct.imp)) {
    if (is.null(pct.reduc)) stop('pct.imp or pct.reduc must be given')
    pct.imp <- 100*(1/(1-pct.reduc/100)-1)
  }
  if (length(control.rate) != length(prop.strat)) 
    stop('length control.rate and prop.strat must be equal')
  a.times <- unique(a.times)
  if (min(a.times)<=0) stop('a.times must be positive')
  if (length(a.times)>1) if (min(diff(a.times)<=1e-7)) 
    stop('a.times must be increasing')
  if (length(a.times)>30) stop('too many analyses')
  if (alpha <= 0 | alpha >= .5) stop('alpha out of range')
  if (alphal >= .5) stop('alphal out of range')
  if (use<1 | use>6) stop('use must be >=1, <=6')
  if (min(c(acc.per,acc.rate,control.rate,pct.imp,prop.strat,p.con))<=0) 
    stop('parameters must be positive')
  if (pct.imp<2) warning('pct.imp must be a % (not a proportion)')
  if (p.con>=1) stop('p.con must be < 1')
  if (use==6 & (oftr<=0 | oftr>=alpha)) stop('Need 0<oftr<alpha')
  if (alphal>0) {
    if (usel<1 | usel>6) stop('usel must be >=1, <=6')
    if (usel==6 & (oftrl<=0 | oftrl>=alphal)) stop('Need 0<oftrl<alphal')
  }
  prop.strat <- prop.strat/sum(prop.strat)
  le <- control.rate/(1+pct.imp/100)
  na <- length(a.times)
  na1 <- na+1
  u2 <- .Fortran('sqopr2',A=as.double(acc.rate),alpha=as.double(alpha),
                 alphal=as.double(alphal),ratimp=as.double(pct.imp),
                 muz=as.double(p.con),pi=as.double(prop.strat),
                 jj=as.integer(length(prop.strat)),lc=as.double(control.rate),
                 le=as.double(le),sa=as.double(acc.per),
                 use=as.integer(use),usel=as.integer(usel),
                 tr=as.double(oftr),trl=as.double(oftrl),
                 kk=as.integer(length(a.times)),rt=as.double(a.times),
                 pt=double(na1),uz=double(na),uzl=double(na),
                 ft0=double(na),ft1=double(na),ac=integer(na),
                 pnu=double(na1),pnl=double(na1),e=double(8),
                 pnu0=double(na1),pnl0=double(na1),eta=double(1),
                 ierr=integer(1),PACKAGE="desmon")[-c(1:15)]
  if (u2$ierr > 0) stop('Boundaries could not be computed')
  w1 <- cbind(u2$pt[-1],u2$ft0,u2$ft1,u2$ac)
  dimnames(w1) <- list(format(round(a.times,2)),c('inf.times','N.fail.H0',
                                                  'N.fail.H1','N.cases'))
  w2 <- cbind(u2$uz,u2$uzl,diff(u2$pnu0),diff(u2$pnl0),diff(u2$pnu),
              diff(u2$pnl))
  if (alphal<=0) w2[,2] <- NA
  dimnames(w2) <- list(dimnames(w1)[[1]],
                       c('Uz','Lz','P.ge.Uz.H0','P.le.Lz.H0','P.ge.Uz.H1',
                         'P.le.Lz.H1'))
  Rej.Probs <- c(size=u2$pnu0[na1],power=u2$pnu[na1],lower.H0=u2$pnl0[na1],
                 lower.H1=u2$pnl[na1])
  w3 <- t(matrix(u2$e,ncol=2))
  dimnames(w3) <- list(c('H0','H1'),c('Real.time','Inf.time',
                                      'N.fail','N.cases'))
  list(call=match.call(),Boundaries=w2,Expected.inf=w1,Rej.Probs=Rej.Probs,
       stop.times=w3,eta=u2$eta)
#   u <- .Fortran('sqopr',A=as.double(acc.rate),alpha=as.double(alpha),
#                 ratimp=as.double(pct.imp),muz=as.double(p.con),
#                 pi=as.double(prop.strat),jj=as.integer(length(prop.strat)),
#                 lc=as.double(control.rate),le=as.double(le),
#                 sa=as.double(acc.per),rt=as.double(a.times),
#                 kk=as.integer(length(a.times)),use=as.integer(use),
#                 pt=double(na1),uz=double(na),ft0=double(na),
#                 ft1=double(na),ac=integer(na),pn0=double(na1),
#                 pn1=double(na1),e=double(8),power=double(1))[-c(1:9,11:12)]
}

#' @export
lr.inf <- function(acc.per,acc.rate,add.fu,control.rate=NULL,pct.imp=NULL,
                   prop.strat=rep(1,length(control.rate)),p.con=.5,
                   control.med=NULL,pct.reduc=NULL) {
  if (is.null(control.rate)) {
    if (is.null(control.med)) stop('control.rate or control.med must be given')
    control.rate <- log(2)/control.med
  }
  if (is.null(pct.imp)) {
    if (is.null(pct.reduc)) stop('pct.imp or pct.reduc must be given')
    pct.imp <- 100*(1/(1-pct.reduc/100)-1)
  }
  if (length(control.rate) != length(prop.strat)) 
    stop('length control.rate and prop.strat must be equal')
  if (min(c(acc.per,acc.rate,control.rate,pct.imp,prop.strat,p.con))<=0) 
    stop('parameters must be positive')
  if (pct.imp<1) warning('pct.imp must be a % (not a proportion)')
  if (p.con>=1) stop('p.con must be < 1')
  prop.strat <- prop.strat/sum(prop.strat)
  theta <- 1 / (1+pct.imp/100)
  le <- control.rate*theta
  d0 <- sum(acc.rate*prop.strat*(acc.per-(exp(control.rate*acc.per)-1)/
                                   exp(control.rate*(acc.per+add.fu))/
                                   control.rate))
  d1 <- sum(acc.rate*prop.strat*(acc.per-(exp(le*acc.per)-1)/
                                   exp(le*(acc.per+add.fu))/le))
  d <- p.con*d0+(1-p.con)*d1
  prop <- p.con*d0/d
  eta <- -sqrt(d*prop*(1-prop))*log(theta)
  c(eta=eta,inf=d)
}

#' @export
seqss <- function(power,acc.per,acc.rate,add.fu,anal.int,first.anal=anal.int,
                  control.rate=NULL,pct.imp=NULL,
                  prop.strat=rep(1,length(control.rate)),
                  alpha=.025,alphal=0,p.con=.5,use=6,
                  usel=6,oftr=alpha/50,oftrl=alphal/50,
                  control.med=NULL,pct.reduc=NULL) {
  if (is.null(control.rate)) {
    if (is.null(control.med)) 
      stop('control.rate or control.med must be given')
    control.rate <- log(2)/control.med
  }
  if (is.null(pct.imp)) {
    if (is.null(pct.reduc)) stop('pct.imp or pct.reduc must be given')
    pct.imp <- 100*(1/(1-pct.reduc/100)-1)
  }
  if (power >= 1 | power <= alpha) stop('power out of range')
  if (use<1 | use>6) stop('use must be >=1, <=6')
  if (min(c(acc.per,acc.rate,control.rate,pct.imp,prop.strat,p.con,anal.int,
            first.anal))<=0) 
    stop('parameters must be positive')
  if (pct.imp<1) warning('pct.imp must be a % (not a proportion)')
  if (acc.per[2]<=acc.per[1]) stop('acc.per must be a positive range')
  FX1 <- function(acc.per) {
    at <- seq(from=first.anal,to=(add.fu+acc.per)*.99,by=anal.int)
    at <- c(at,add.fu+acc.per)
    seqopr(acc.per,acc.rate,control.rate,pct.imp,at,prop.strat,alpha,alphal,
           p.con,use,usel,oftr,oftrl)$Rej.Probs[2]-power
  }
  z <- uniroot(FX1,acc.per)
  at <- seq(from=first.anal,to=(add.fu+z$root)*.99,by=anal.int)
  at <- c(at,add.fu+z$root)
  u <- seqopr(z$root,acc.rate,control.rate,pct.imp,at,prop.strat,alpha,
              alphal,p.con,use,usel,oftr,oftrl)
  u$acc.per <- z$root
  u
}
