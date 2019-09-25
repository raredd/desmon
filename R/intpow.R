#' Approximate Power of Cox Model Interaction Tests
#' 
#' Approximates the power of the large sample partial likelihood tests for a
#' treatment by biologic marker interaction test in a Cox model using the
#' normal approximation, assuming exponential distributions.
#' 
#' @aliases intpow hazsf expinf siminf
#' 
#' @details
#' For a Cox model with two binary factors (x=treatment and z=biologic marker),
#' \code{intpow} computes the power of the large sample partial likelihood
#' tests for no interaction. Constant hazards are assumed within each of the 4
#' covariate combinations.
#' 
#' Let \code{p(x,z)} be the sample proportion and \code{h(x,z)} be the hazard 
#' for treatment \code{x} in marker group \code{z}. Treatment is assumed 
#' independent of the marker, so the proportion of the sample in the 4 
#' combinations is \code{p(0,0)=(1-prx)*(1-prz),
#' p(0,1)=(1-prx)*prz, p(1,0)=prx*(1-prz), p(1,1)=prx*prz}.
#' 
#' The hazards can either be given explicitly (argument \code{hazs}) or 
#' through the ratios and the overall average hazard. The ratios are defined by
#' \code{ri=h(1,1)*h(0,0)/(h(0,1)*h(1,0))},
#' \code{rt=((1-prz)*h(1,0)+prz*h(1,1))/((1-prz)*h(0,0)+prz*h(0,1))}, and
#' \code{rf=((1-prx)*h(0,1)+prx*h(1,1))/((1-prx)*h(0,0)+prx*h(1,0))}.
#' 
#' Given a set of ratios and an average hazard \code{l0}, \code{hazsf} computes
#' the hazards in the individual groups to give the specified ratios and to 
#' satisfy the constraint
#' \code{p(0,0)*h(0,0)+p(0,1)*h(0,1)+p(1,0)*h(1,0)+p(1,1)*h(1,1)=l0}.
#' \code{hazsf} is called by \code{intpow} if \code{hazs} is not specified.
#' 
#' The calculations assume a clinical trial setting with accrual uniform over
#' the period \code{(0,acc.per)} and an additional \code{add.fu} units of
#' follow-up after completion of accrual. Thus censoring is assumed uniform on
#' \code{(add.fu, add.fu+acc.per)}.
#' 
#' \code{expinf} computes the expected information matrix (per subject) from 
#' the Cox partial likelihood containing the covariates \code{x}, \code{z}, 
#' and \code{x*z}.
#' 
#' \code{siminf} calculates the Wald statistic for the interaction test and 
#' the partial likelihood variance estimate of the variance of the interaction 
#' coefficient for a single simulated sample from this model. 
#' 
#' \code{intpow} calculates the power for the Wald test for the null hypothesis
#' that coefficient of the interaction term = 0 using the large sample normal 
#' approximation. In large samples, this test is equivalent to the score and 
#' partial likelihood ratio tests for this hypothesis.
#' 
#' @usage
#' intpow(n, prx, prz, ri, rt, rf, l0, acc.per, add.fu, 
#'        alpha2 = .05, hazs = NULL)
#' hazsf(prx, prz, ri, rt, rf, l0)
#' expinf(prx, prz, hazs, acc.per, add.fu)
#' siminf(prx, prz, hazs, acc.per, add.fu, n)
#' 
#' @param n Sample size
#' @param prx Proportion randomized to treatment 1 (x=1)
#' @param prz Proportion positive for the marker (z=1)
#' @param ri Interaction ratio: the treatment hazard ratio (treatment 1 /
#' treatment 0) in the marker positive group divided by the treatment hazard
#' ratio in the marker negative group
#' @param rt Average treatment hazard ratio
#' @param rf Average marker hazard ratio (positive/negative)
#' @param l0 Overall average hazard
#' @param acc.per Number of time units of accrual
#' @param add.fu Number of time units of follow-up after the end of accrual
#' @param alpha2 Two-sided type I error rate for the interaction test
#' @param hazs Vector giving the constant hazard rates in the 4 groups, in the
#' order \code{(x,z)=(0,0),(0,1),(1,0),(1,1)}
#' @return \code{intpow} returns a list with components: \item{power }{ The
#' power of the interaction test} \item{hazards }{ The vector of hazards from
#' \code{hazsf}} \item{var}{ The inverse expected information matrix (per
#' sample)} \item{nevents}{ The expected number of events in each of the
#' \code{x} by \code{z} combinations}
#' 
#' \code{hazsf} returns the vector of hazards meeting the constraints on the
#' ratios and average hazard, in the order \code{(x,z)=(0,0),(0,1),(1,0),(1,1)}
#' 
#' \code{expinf} returns the expected information matrix per subject
#' 
#' \code{siminf} returns a vector containing the Wald test statistic and the
#' estimated variance of the interaction term
#' @keywords design survival
#' @examples
#' intpow(n = 2800, prx = .5, prz = .25, ri = .6, rt = .8, rf = 2, 
#'        l0 = .045, acc.per = 2.5, add.fu = 3)
#' 
#' @export intpow

#' @export
intpow <- function(n, prx, prz, ri, rt, rf, l0, acc.per, add.fu, 
                   alpha2 = .05, hazs = NULL) {
  # power for treatment by biologic marker interaction tests
  # n=total sample size
  # prx= proportion randomized to first treatment
  # prz= proportion marker positive
  # ri=interaction treatment hazard ratio ratio: 
  # [(x=1,z=1)/(x=0,z=1)]/[(x=1,z=0](x=0,z=0)]
  # rt=average treatment hazard ratio (1/0)
  # rf=average marker hr (1/0)
  # l0=average hazard
  # acc.per=time in years to accrue n patients
  # add.fu= additional years follow-up after end of accrual
  # alpha2=two-sided type I error
  if (is.null(hazs)) hazs <- hazsf(prx=prx,prz=prz,ri=ri,rt=rt,rf=rf,l0=l0)
  else ri <- hazs[4]*hazs[1]/(hazs[2]*hazs[3])
  u1 <- expinf(prx=prx,prz=prz,hazs=hazs,acc.per=acc.per,add.fu=add.fu)
  u1i <- solve(u1)/n
  pxz <- c((1-prx)*c((1-prz),prz),prx*c((1-prz),prz))
  nev <- rep(0,length(pxz))
  for (i in 1:length(pxz)) nev[i] <- nevents(n*pxz[i]/acc.per,acc.per,
                                             add.fu,haz=hazs[i])[1]
  names(nev) <- names(hazs)
  list(power=1-pnorm(qnorm(1-alpha2/2)-abs(log(ri))/sqrt(u1i[3,3])),
       hazards=c(hazs),var=u1i,nevents=nev)
}

#' @export
siminf <- function(prx,prz,hazs,acc.per,add.fu,n) {
  pxz <- round(c((1-prx)*c((1-prz),prz),prx*c((1-prz),prz))*n)
  grp <- rep(1:4,pxz)
  cov <- (rbind(c(0,0,0),c(0,1,0),c(1,0,0),c(1,1,1)))[grp,]
  ft <- rexp(n)/hazs[grp]
  ct <- runif(n)*(acc.per)+add.fu
  status <- ifelse(ct<ft,0,1)
  ft <- pmin(ft,ct)
  z <- coxph(Surv(ft,status)~cov[,1]+cov[,2]+cov[,3])
  c(abs(z$coef)[3]/sqrt(z$var[3,3]),z$var[3,3])
}

#' @export
expinf <- function(prx,prz,hazs,acc.per,add.fu) {
  # Cox model expected information per subject
  # constant hazards in each group
  pxz <- c((1-prx)*c((1-prz),prz),prx*c((1-prz),prz))
  u2 <- u <- matrix(0,3,3)
  g <- function(t,ap,add.fu) {#censoring survivor function
    pmin(1,pmax((ap+add.fu-t)/ap,0))
  }
  f1 <- function(x,hazs,pxz,comp,ap,add.fu,g) {
    pxz[comp]*hazs[comp]*exp(-x*hazs[comp])*g(x,ap,add.fu)
  }
  f2 <- function(x,hazs,pxz,c1,c2,ap,add.fu,g) {
    w <- hazs*pxz
    u <- cbind(w[1]*exp(-x*hazs[1]),w[2]*exp(-x*hazs[2]),w[3]*exp(-x*hazs[3]),
               w[4]*exp(-x*hazs[4]))*g(x,ap,add.fu)
    u2 <- c(u %*% rep(1,4))
    u3 <- if (length(c1)>1) c(u[,c1] %*% rep(1,length(c1))) else c(u[,c1])
    u4 <- if (length(c2)>1) c(u[,c2] %*% rep(1,length(c2))) else c(u[,c2])
    ifelse(u2>0,u3*u4/u2,0)
  }
  u[1,1] <- u[2,2] <- u[1,2] <- u[1,3] <- u[2,3] <- u[3,3] <- 
    do.call('integrate',list(f=f1,lower=0,upper=acc.per+add.fu,hazs=hazs,
                             pxz=pxz,comp=4,ap=acc.per,add.fu=add.fu,g=g))[[1]]
  u[1,1] <- u[1,1]+
    do.call('integrate',list(f=f1,lower=0,upper=acc.per+add.fu,hazs=hazs,
                             pxz=pxz,comp=3,ap=acc.per,add.fu=add.fu,g=g))[[1]]
  u[2,2] <- u[2,2]+
    do.call('integrate',list(f=f1,lower=0,upper=acc.per+add.fu,hazs=hazs,
                             pxz=pxz,comp=2,ap=acc.per,add.fu=add.fu,g=g))[[1]]
  #  print(c(failprob=(u[1,1]+u[2,2]-u[3,3]+
  #                   do.call('integrate',
  #                           list(f=f1,lower=0,upper=acc.per+add.fu,hazs=hazs,
  #                                pxz=pxz,comp=1,ap=acc.per,
  #                                add.fu=add.fu,g=g))[[1]])))
  u2[1,1] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(3,4),c2=c(3,4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u2[1,2] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(3,4),c2=c(2,4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u2[1,3] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(3,4),c2=c(4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u2[2,2] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(2,4),c2=c(2,4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u2[2,3] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(2,4),c2=c(4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u2[3,3] <- do.call('integrate',
                     list(f=f2,lower=0,upper=acc.per+add.fu,hazs=hazs,
                          pxz=pxz,c1=c(4),c2=c(4),ap=acc.per,
                          add.fu=add.fu,g=g))[[1]]
  u <- u-u2
  u[2:3,1] <- u[1,2:3]
  u[3,2] <- u[2,3]
  u
}

#' @export
hazsf <- function(prx,prz,ri,rt,rf,l0) {
  # calculate hazards for 4 groups such that the interaction ri is as 
  # specified and the marginal ratios and average rate constraints are met.
  # called by intpow
  pxz <- c((1-prx)*c((1-prz),prz),prx*c((1-prz),prz))
  u <- rbind(c((1-prx)*c(-rf,1),prx*c(-rf,1)),
             c(-rt*c(1-prz,prz),1-prz,prz),pxz)
  rhs <- c(0,0,l0)
  a1 <- svd(u,3,4)
  ns <- a1$v[,4]
  ps <- a1$v[,1:3] %*% ((t(a1$u) %*% rhs)/a1$d)
  aa <- ns[1]*ns[4]-ri*ns[2]*ns[3]
  cc <- ps[1]*ps[4]-ri*ps[2]*ps[3]
  bb <- ps[1]*ns[4]+ps[4]*ns[1]-ri*(ns[2]*ps[3]+ns[3]*ps[2])
  u2 <- (-bb+c(-1,1)*sqrt(bb*bb-4*aa*cc))/(2*aa)
  #  print(u2)
  h1 <- ps+u2[1]*ns
  if (min(h1) <= 0) h1 <- ps+u2[2]*ns
  #  u3 <- h1*pxz
  #  print(c(sum(u3),sum(u3[3:4])/sum(u3[1:2]),sum(u3[c(2,4)])/sum(u3[c(1,3)]),
  #        u3[1]*u3[4]/(u3[2]*u3[3]))) # checks constraints
  names(h1)=c('x0z0','x0z1','x1z0','x1z1')
  h1
}