c Version 1.0                     seqopr2                  
c
c		Original Author: Kyungmann Kim, Ph.D.
c		Language: Fortran 77
c
c ======================== Program Description ========================
c
c	Operating Characteristic for Group Sequential Clinical Trials
c	with Censored Survival Data, Adjusting for Stratification and
c	with Possibly Unequal Patient Allocation Between Control and Rx
c
c for a group sequential clinical trial with censored survival data
c adjusting for stratification and possibly unequal patient allocation,
c this program computes the expected stopping time in the information
c time, the real chronological time, the number of failures, and the
c number of patients.
c this program can be used to investigate the selected group sequential
c design for its operating characteristics.
c
c See Kim and Cheuvart (1990).  Study Duration for Group Sequential
c Clinical Trials with Censored Survival Data Adjusting for Stratification.
c
c
c =====================================================================

      subroutine sqopr2(A,alpha,alphal,ratimp,muz,pi,jj,lc,le,sa,
     $     use,usel,tr,trl,kk,rt,pt,uz,uzl,ft0,ft1,ac,pnu,pnl,
     $     e,pnu0,pnl0,eta,ierr)

c *********************************************************************
c compute the operating characteristics of group sequential clinical
c trials with censored survival data, adjusting for stratification and
c with possibly unequal patient allocations between treatments
c *********************************************************************
c ====================== Key Design Parameters ========================
c
c      kk := maximum number of repeated analyss
c     use := indicators for use function
c    usel := indicators for use function for rci lower boundary
c      rt := real chronological times of repeated analyses
c   alpha := significance level or type I error 
c  alphal := one-sided rci lower boundary alpha
c  ratimp := percent improvement in survival from control to treatment
c     muz := fraction of patients randomized to control
c      jj := number of strata
c      pi := fraction of patients in strata
c       A := accrual rate per time unit
c      lc := control hazard rates (1 to # strata)
c      le := experimental hazard rates (1 to # strata)
c      sa := accrual duration in time units
c  (note, must have le=lc/(1+ratimp/100)
c
c   pt = information times at each analysis
c   uz = upper boundary at each analysis
c   ft0 = expected # failures under H0 at each analysis
c   ft1 = expected # failures under H1 at each analysis
c   ac = # patients accrued at each analysis
c   pnl = cumulative lower stopping probability at each analysis H1
c   pnu = cumulative upper stopping probability at each analysis H1
c   pnl0 = cumulative lower stopping probability at each analysis H0
c   pnu0 = cumulative upper stopping probability at each analysis H0
c  ept0 = exp inf time at termination under H0
c  ert0 = exp real time at term under H0
c  eft0 = exp # failures at term under H0
c  eac0 = exp # patients at term under H0
c  ept1 ert1 eft1 eac1 = same quantities under H1
c  put in array e
      integer jj,kk,use,usel,ierr
      real*8 A,alpha,alphal,ratimp,muz,tr,trl
      real*8 pi(jj),lc(jj),le(jj)
      real*8 sa,rt(kk),e(8)

c yn, nyn upper and lower boundaries on partial sum scale
      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

      integer k,ac(kk)
      real*8 eta,theta
      real*8 Emax,prop,pnu0(0:kk),pnl0(0:kk)
      real*8 pnu(0:kk),pnl(0:kk),uz(kk),uzl(kk)
      real*8 pt(0:kk),ept0,ept1,ft0(kk),ft1(kk),eft0,eft1
      real*8 ert0,ert1,eac0,eac1,tmp,tmp2
      real*8 nfa,nff

      h = 5.d-2
      pt(0) = 0.d0
      pnu(0) = 0.d0
      pnl(0) = 0.d0
      pnu0(0) = 0.d0
      pnl0(0) = 0.d0

c theta= hazard ratio exp/control (<1)
      theta = 1.d0 / (1.d0+ratimp/100.d0)

      do 20 k=1,kk
         if(rt(k) .lt. sa) then
            ft0(k) = nfa(rt(k),lc,A,jj,pi)
            ft1(k) = muz*ft0(k)+(1-muz)*nfa(rt(k),le,A,jj,pi)
         else
            ft0(k) = nff(rt(k),lc,A,sa,jj,pi)
            ft1(k) = muz*ft0(k)+(1-muz)*nff(rt(k),le,A,sa,jj,pi)
         endif
 20   continue

c maximum number of expected failures under the alternative hypothesis

      Emax = ft1(kk)
c      Emax = muz*nff(rt(kk),lc,A,sa,jj,pi)+
c     +     (1.d0-muz)*nff(rt(kk),le,A,sa,jj,pi)

c determine "information times" based on expected number of failures
c under "null" hypothesis -- now using alt hypothesis

      do 21 k=1,kk
         pt(k)=ft1(k)/ft1(kk)
c         pt(k) = ft0(k) / ft0(kk)
c         ft1(k) = pt(k) * Emax
 21   continue

c determine accrual time based on information time

      do 22 k=1,kk
         if(rt(k) .lt. sa) then
            ac(k) = A * rt(k) + 0.5d0
         else
            ac(k) = A * sa + 0.5d0
         endif
 22   continue

c compute drift parameter; prop=expected proportion of failures on control arm under Ha

      prop = muz*ft0(kk) / Emax
      eta = -dsqrt(Emax*prop*(1.d0-prop))*dlog(theta)

c generate one-sided group sequential boundary using use function for rci

      do 10 k=1,kk
         nyn(k)=-5
 10   continue 
      if (alphal.gt.0) then
         call ones(kk,alphal,pt,pnl,usel,trl,ierr)
         if (ierr.eq.1) return
c compute approx rci lower boundary
         do 810 k=1,kk
            nyn(k)=-yn(k)+eta*pt(k)
 810     continue 
      endif

c compute upper boundary

      call ones(kk,alpha,pt,pnl,use,tr,ierr)
      if (ierr.eq.1) return

c standardize boundaries

      do 30 k=1,kk
         tmp = dsqrt(pt(k))
         if (nyn(k).gt.yn(k)) nyn(k)=yn(k)
         uz(k) = yn(k)/tmp
         uzl(k) = nyn(k)/tmp
 30   continue

c boundary crossing probs under null

      call pete2(kk,pt,pnu0,pnl0,0.d0,2,tmp,tmp2,ierr)
      if (ierr.eq.1) return
c compute power 

      call pete2(kk,pt,pnu,pnl,eta,2,tmp,tmp2,ierr)
      if (ierr.eq.1) return

      if (alphal.le.0) then
         do 33 k=1,kk
            pnl0(k)=0
            pnl(k)=0
 33      continue 
      endif

c compute expected stopping in inf time, real time, failures, and accruals

      ept0 = 0.d0
      ept1 = 0.d0
      eac0 = 0.d0
      eac1 = 0.d0
      ert0 = 0.d0
      ert1 = 0.d0
      eft0 = 0.d0
      eft1 = 0.d0
      do 40 k=1,kk-1
         tmp=pnu0(k)-pnu0(k-1)+pnl0(k)-pnl0(k-1)
         ept0 = ept0 + tmp * pt(k)
         eac0 = eac0 + tmp * ac(k)
         ert0 = ert0 + tmp * rt(k)
         eft0 = eft0 + tmp * ft0(k)
         tmp=pnu(k)-pnu(k-1)+pnl(k)-pnl(k-1)
         ept1 = ept1 + tmp * pt(k)
         eac1 = eac1 + tmp * ac(k)
         ert1 = ert1 + tmp * rt(k)
         eft1 = eft1 + tmp * ft1(k)
 40   continue
      tmp=1-pnu0(kk-1)-pnl0(kk-1)
      ept0 = ept0 + tmp * pt(kk)
      eac0 = eac0 + tmp * ac(kk)
      ert0 = ert0 + tmp * rt(kk)
      eft0 = eft0 + tmp * ft0(kk)
      tmp=1-pnu(kk-1)-pnl(kk-1)
      ept1 = ept1 + tmp * pt(kk)
      eac1 = eac1 + tmp * ac(kk)
      ert1 = ert1 + tmp * rt(kk)
      eft1 = eft1 + tmp * ft1(kk)
      e(1)=ert0
      e(2)=ept0
      e(3)=eft0
      e(4)=eac0
      e(5)=ert1
      e(6)=ept1
      e(7)=eft1
      e(8)=eac1

      return
      end

      subroutine sqopr3(kk,pt,eta,uz,uzl,pnu,pnl,ierr)
c
c this version just computes the probability of crossing a specified boundary
c
c      kk := maximum number of repeated analyses
c   pt = information times at each analysis
c   eta = brouwnian motion drift parameter (related to alternative)
c   uz = upper boundary at each analysis
c   uzl = lower boundary at each analysis
c   pnl = cumulative lower stopping probability at each analysis
c   pnu = cumulative upper stopping probability at each analysis

      integer kk,ierr,k
      real*8 eta
      real*8 pnu(0:kk),pnl(0:kk),uz(kk),uzl(kk)
      real*8 pt(0:kk),tmp,tmp2

c yn, nyn upper and lower boundaries on partial sum scale
      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

      h = 5.d-2
      pt(0) = 0.d0
      pnu(0) = 0.d0
      pnl(0) = 0.d0

c convert boundaries to partial sum scale

      do 30 k=1,kk
         tmp = dsqrt(pt(k))
         yn(k)=uz(k)*tmp
         nyn(k)=uzl(k)*tmp
         if (nyn(k).gt.yn(k)) nyn(k)=yn(k)
c         uz(k) = yn(k)/tmp
c         uzl(k) = nyn(k)/tmp
 30   continue

c compute boundary crossing probs

      call pete2(kk,pt,pnu,pnl,eta,2,tmp,tmp2,ierr)
c      if (ierr.eq.1) return

      return
      end

      subroutine pete2(k,pt,pn,pnl,eta,side,pete,petl,ierr)

c This program computes the probability of crossing upper or either
c boundaries using the iterative numerical integration method suggested
c by Armitage et al. (1969) and McPherson and Armitage (1971).

c a calling sequence

      integer k,side,ierr
      real*8 pt(0:k),pn(0:k),pnl(0:k),eta,pete,petl

c local definitions

      integer i,j,lnp2
      real*8 fnn,dnni,sumn,tempn,qneg
      real*8 fnp,dnpi,sump,tempp,qpos
      real*8 eps,mean,varian
      real*8 fcab,dgaudf

c common definitions

      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

      eps = 1.d-7
      pete = 0.d0
      petl = 0.d0

      do 600 i=1,k
         varian = pt(i) - pt(i-1)
         mean = eta * varian

c compute fcab for later use in iterative integration
         ln(i) = (yn(i) - nyn(i)) / h
         if(i .lt. k) then
            lnp2 = ln(i) + 2
            if (lnp2.gt.500) then
               ierr=1
               return
            endif
            dnpi = nyn(i)
            do 300 j=1,lnp2
               if(j .eq. lnp2) dnpi = yn(i)
               fn(j,i) = fcab(dnpi,i,mean,varian)
               dnpi = dnpi + h
 300        continue 
         endif

         if(i .eq. 1) then
            qpos = 1.d0 - dgaudf((yn(i)-mean)/dsqrt(varian))
            if(side .gt. 1) qneg = dgaudf((nyn(i)-mean)/dsqrt(varian))
            goto 575
         endif

c integrate fcab from yn(i) to infinity for crossing upper boundary

         dnpi = yn(i)
         sump = fcab(dnpi,i,mean,varian)

 400     dnpi = dnpi + h
         fnp = fcab(dnpi,i,mean,varian)
         sump = sump + 2.d0 * fnp

c increase counter and grid points
            
         tempp = h * fnp
         if (tempp.gt.eps) goto 400

c far enough out to +/- infinity

         sump = sump - fnp
         qpos = h * sump / 2.d0
         
         if (side .eq. 1) goto 580

c integrate fcab from nyn(i) to -infinity for crossing lower boundary

 450     dnni = nyn(i)
         sumn = fcab(dnni,i,mean,varian)

 500     dnni = dnni - h
         fnn = fcab(dnni,i,mean,varian)
         sumn = sumn + 2.d0 * fnn

c increase counter and grid points

         tempn = h * fnn
         if (tempn.gt.eps) goto 500

c far enough out to +/- infinity
            
         sumn = sumn - fnn
         qneg = h * sumn / 2.d0

c accumulate tail probabilities

 575     continue
c         if(side .gt. 1) qpos = qpos + qneg

 580     pete = pete + qpos
         pn(i) = pete
         petl = petl + qneg
         pnl(i) = petl
 600  continue 
      
      return
      end
