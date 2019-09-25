#' Computes the power of the logrank test
#' 
#' Computes the power of the two-group logrank test for arbitrary failure time
#' distributions in a standard clinical trials setting.
#' 
#' @details
#' The calculations assume that \code{n=acc.rate*acc.per} patients are entered
#' uniformly over the period \code{[0,acc.per]}, with follow-up continuing for
#' an addition \code{add.fu} time units, so censoring will be uniform on
#' \code{[add.fu,add.fu+acc.per]}.  The failure probability under the specified
#' failure distribution is then computed, as is the expected number of
#' failures, \code{n*fail.prob}.
#' 
#' \code{hazcon} must be the hazard function corresponding to the survivor
#' function in \code{survcon} (and similarly for \code{haztst} and
#' \code{survtst}).  The program does not check for consistency of these
#' functions, though.  One way to check would be to compare \code{survcon(x)}
#' to \code{exp(-integrate(hazcon,0,x)$value)} for various values \code{x}.
#' The default functions are for the exponential cure rate model, which reduces
#' to the two-sample exponential when the cure fractions are 0.
#' 
#' If custom \code{hazcon}, \code{survcon}, \code{haztst}, and \code{survtst}
#' functions are given, the arguments should match those in the default
#' functions (ie \code{x,lc,pc,...} for \code{hazcon}, etc).  The values of
#' \code{lc}, \code{pc}, \code{lt}, and \code{pt} are specified in the
#' \code{control.rate}, \code{control.cure}, \code{test.rate}, and
#' \code{test.cure} arguments and passed through to these functions.
#' 
#' The calculations are performed by using the \code{integrate} function to
#' approximate the expectations of the logrank score and the logrank variance
#' estimator and the variance of the logrank score under the specified
#' conditions, and then using a normal approximation to the distribution of the
#' logrank statistic.
#' 
#' This function only computes power for a single analysis, and does not
#' consider group sequential tests.
#' 
#' Caution: consistent time units must be used for all quantities; eg if the
#' accrual rate is given in patients/month, then the hazard rates must be in
#' units of events/month and \code{acc.per} and \code{add.fu} must be in
#' months.
#' 
#' @usage
#' powlgrnk(acc.per, acc.rate, add.fu, alpha = .025, p.con = .5,
#'          hazcon = function(x,lc,pc,...) 
#'                      {u <- (1 - pc) * exp(-lc * x); lc * u / (pc + u)},
#'          survcon = function(x, lc, pc, ...) pc+(1-pc)*exp(-lc*x),
#'          haztst = function(x, lt, pt, ...) 
#'                      {u <- (1 - pt) * exp(-lt * x); lt * u / (pt + u)},
#'          survtst = function(x, lt, pt, ...) pt+(1-pt)*exp(-lt*x),
#'          control.rate = NULL, test.rate = NULL, control.cure = 0,
#'          test.cure = 0, ...)
#' 
#' @param acc.per Planned duration of accrual
#' @param acc.rate Number of patients expected to be entered per time unit
#' @param add.fu Additional follow-up between the end of accrual and the time
#' of analysis
#' @param alpha The one-sided type I error of the test
#' @param p.con The proportion randomized to the control arm
#' @param hazcon A function evaluating the control group hazard function at a
#' vector of times
#' @param survcon A function evaluating the control group survivor function at
#' a vector of times
#' @param haztst A function evaluating the experimental or test group hazard
#' function at a vector of times
#' @param survtst A function evaluating the experimental or test group survivor
#' function at a vector of times
#' @param control.rate Exponential hazard rate in control group for non-cured
#' @param test.rate Exponential hazard rate in test group for non-cured
#' @param control.cure Cure fraction in control group
#' @param test.cure Cure fraction in test group
#' @param ... additional arguments passed to survival and hazard functions
#' 
#' @return Returns a vector giving the power of the test under the specified
#' conditions (\code{power}), the total sample size (\code{n}) and the expected
#' number of events (\code{nd}).
#' 
#' @seealso
#' \code{\link{powlgrnk6}}
#' 
#' @keywords survival design
#' 
#' @examples
#' # Exponential distributions
#' powlgrnk(5,200,3,control.rate=.1,test.rate=.075)
#' 
#' #       power            n           nd 
#' #   0.7925548 1000.0000000  375.5712999 
#' 
#' # Cure rate
#' powlgrnk(3,200,3,control.rate=log(2)/3,test.rate=log(2)/4,
#'   control.cure=.3,test.cure=.4)
#'   
#' #      power           n          nd 
#' #  0.8962665 600.0000000 230.7956591
#' 
#' # Exponential cure-rate with proportional hazards alternative
#' ht <- function(x,rat,...) rat*hc(x,...)
#' st <- function(x,rat,...) sc(x,...)^rat
#' sc <- function(x,pic,lic,...) pic+(1-pic)*exp(-lic*x)
#' hc <- function(x,pic,lic,...) {u <- (1-pic)*exp(-lic*x); lic*u/(pic+u)}
#' powlgrnk(5,200,3,hazcon=hc,survcon=sc,haztst=ht,survtst=st,lic=log(2)/3,
#'   pic=.3,rat=.5)
#'   
#' #       power            n           nd 
#' #   0.8564817 1000.0000000  446.0797311 
#' 
#' @export

powlgrnk <- function(acc.per, acc.rate, add.fu, alpha = .025, p.con = .5,
                     hazcon = function(x,lc,pc,...) 
                       {u <- (1 - pc) * exp(-lc * x); lc * u / (pc + u)},
                     survcon = function(x, lc, pc, ...) pc+(1-pc)*exp(-lc*x),
                     haztst = function(x, lt, pt, ...) 
                       {u <- (1 - pt) * exp(-lt * x); lt * u / (pt + u)},
                     survtst = function(x, lt, pt, ...) pt+(1-pt)*exp(-lt*x),
                     control.rate = NULL, test.rate = NULL, control.cure = 0,
                     test.cure = 0, ...) {
    
  ### Computes power and expected #events (nd) for logrank using large sample
  ###   mean and variance formulas 
  ### survcon and hazcon are functions evaluating the the control group 
  ###   survival and hazard functions
  ### survtst and haztst similarly for the test group
  ### acc.per=accrual period
  ### add.fu=additional follow-up after completion of accrual
  ### alpha=one-sided type I error,
  ### p.con=proportion randomized to control
  risk <- function(t,c.rng,surv,...) ifelse (t<=c.rng[1],1,ifelse(
     t>=c.rng[2],0,(c.rng[2]-t)/((c.rng[2]-c.rng[1]))))*surv(t,...) 
     #P(still at risk at t)
  mn <- function(t,p.con,c.rng,risk,survcon,hazcon,survtst,haztst,
      lc,lt,pc,pt,...) {#mean of logrank score
    y1 <- p.con*risk(t,c.rng,survcon,lc=lc,pc=pc,...)
    y2 <- (1-p.con)*risk(t,c.rng,survtst,lt=lt,pt=pt,...)
    y <- y1+y2
    ifelse(y>0,y1*y2*(hazcon(t,lc=lc,pc=pc,...)-haztst(t,lt=lt,pt=pt,...))/y,0)
  }
  vs <- function(t,p.con,c.rng,risk,survcon,hazcon,survtst,haztst,
      lc,lt,pc,pt,...) {#expected value of null variance estimator
    y1 <- p.con*risk(t,c.rng,survcon,lc=lc,pc=pc,...)
    y2 <- (1-p.con)*risk(t,c.rng,survtst,lt=lt,pt=pt,...)
    y <- y1+y2
    ifelse(y>0,(y1*y2/y)^2*(hazcon(t,lc=lc,pc=pc,...)/y2+haztst(t,lt=lt,pt=pt,...)/y1),0)
  }
  sig2 <- function(t,p.con,c.rng,risk,survcon,hazcon,survtst,haztst,
      lc,lt,pc,pt,...) {#variance of logrank score
    y1 <- p.con*risk(t,c.rng,survcon,lc=lc,pc=pc,...)
    y2 <- (1-p.con)*risk(t,c.rng,survtst,lt=lt,pt=pt,...)
    y <- y1+y2
    ifelse(y>0,(y1*y2/y)^2*(haztst(t,lt=lt,pt=pt,...)/y2+hazcon(t,lc=lc,pc=pc,...)/y1),0)
  }
  pf <- function(t,c.rng,risk,surv,haz,...) haz(t,...)*risk(t,c.rng,surv,...)
  c.rng <- c(add.fu,acc.per+add.fu)
  r1 <- integrate(mn,0,c.rng[2],p.con=p.con,c.rng=c.rng,risk=risk,
    survcon=survcon,survtst=survtst,hazcon=hazcon,haztst=haztst,
    lc=control.rate,lt=test.rate,pc=control.cure,pt=test.cure,...)[[1]]
  r2 <- integrate(vs,0,c.rng[2],p.con=p.con,c.rng=c.rng,risk=risk,
    survcon=survcon,survtst=survtst,hazcon=hazcon,haztst=haztst,
    lc=control.rate,lt=test.rate,pc=control.cure,pt=test.cure,...)[[1]]
  r3 <- integrate(sig2,0,c.rng[2],p.con=p.con,c.rng=c.rng,risk=risk,
    survcon=survcon,survtst=survtst,hazcon=hazcon,haztst=haztst,
    lc=control.rate,lt=test.rate,pc=control.cure,pt=test.cure,...)[[1]]
  r4 <- integrate(pf,0,c.rng[2],c.rng=c.rng,risk=risk,surv=survcon,
    haz=hazcon,lc=control.rate,pc=control.cure,...)[[1]]
  r5 <- integrate(pf,0,c.rng[2],c.rng=c.rng,risk=risk,surv=survtst,
    haz=haztst,lt=test.rate,pt=test.cure,...)[[1]]
  crit <- qnorm(1-alpha)
  n <- acc.rate * acc.per
  power <- 1-pnorm(crit*sqrt(r2/r3)-r1*sqrt(n/r3))
  c(power=power,n=n,nd=n*(p.con*r4+(1-p.con)*r5))
}
