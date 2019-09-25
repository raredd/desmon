      subroutine ones (k,alpha,pt,pn,star,tr,ierr)

c This subroutine generates one-sided boundaries for group sequential
c tests of H0:eta<=0 vs H1:eta>0 according to the Type I Error
c Spending Rate Function alpha * in Lan and DeMets (1983) by using
c the repeated numerical integration method suggested by Armitage
c et al. (1969).  

c    a calling sequence

      integer k,star,ierr
      real*8 alpha,pt(0:k),pn(0:k),tr

c    common definitions

      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

c    local definitions

      integer i,j,lnp2,flag
      real*8 a,b,dnp,sum,temp,fold,fnew
      real*8 pd,eps,mean,var
      real*8 fcab,falsi3,dgauin
      real*8 alpha1,alpha2,alpha3,alpha4,alpha5
      
      ierr=0
      if (star.gt.5) then 
	 call ones6(k,alpha,pt,pn,tr,ierr)
	 return
      endif
      eps = 1.d-7
      mean = 0.d0

      do 300 i=1,k

c       probability of crossing upper boundary by pt(i)

         if(star .eq. 1) then
            pn(i) = alpha1(pt(i),alpha)
         elseif(star .eq. 2) then
            pn(i) = alpha2(pt(i),alpha)
         elseif(star .eq. 3) then
            pn(i) = alpha3(pt(i),alpha)
         elseif(star .eq. 4) then
            pn(i) = alpha4(pt(i),alpha)
         else
            pn(i) = alpha5(pt(i),alpha)
         endif

         pd = pn(i) - pn(i-1)
         var = pt(i) - pt(i-1)

         if(pd .le. 0.d0) then
            yn(i) = 5.d0
c          nyn(i) = -5.d0
            goto 175
         endif

         if(i .eq. 1) then
            yn(i) = dgauin(1.d0-pd)*dsqrt(var) + mean
            goto 175
         endif

c       integrate fcab from infinity to yn(i)

         sum = 0.d0
         dnp = 5.d0
         fold = fcab(dnp,i,mean,var)
 100     dnp = dnp - h
         fnew = fcab(dnp,i,mean,var)
         temp = h * (fold+fnew) / 2.d0
         sum = sum + temp

c       check if upper limit is large enough

         if(sum .ge. pd) goto 150
         if (dnp.lt.0) then
            ierr=1
            return
         endif
         fold = fnew
         goto 100

c       initialize for falsi3

 150     dnp = dnp + h
         sum = sum - temp
         a = 0.d0
         b = h
         flag = -1

c call falsi3 for interpolation

         yn(i) = falsi3(flag,dnp,fold,i,a,b,sum,temp,pd,mean,var)

c determine the number of grids from -infinity to yn(i)

c175     nyn(i) = -5.d0
 175     ln(i) = (yn(i) - nyn(i)) / h

         if(i .eq. k) goto 400

c       compute fcab for later use

         lnp2 = ln(i) + 2
         if (lnp2.gt.500) then
            ierr=1
            return
         endif
         dnp = nyn(i)
         do 200 j=1,lnp2
            if(j .eq. lnp2) dnp = yn(i)
            fn(j,i) = fcab(dnp,i,mean,var)
            dnp = dnp + h
 200     continue
 300  continue

 400  return
      end

      subroutine ones6 (k,alpha,pt,pn,tr,ierr)

c For truncated O-F

c    a calling sequence

      integer k,ierr
      real*8 alpha,pt(0:30),pn(0:30),tr

c    common definitions

      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

c    local definitions

      integer i,j,lnp2,flag,flag2
      real*8 a,b,dnp,sum,temp,fold,fnew
      real*8 pd,eps,mean,var
      real*8 fcab,falsi3,dgauin
      real*8 alpha1
      real*8 qs
c qs=truncated boundary on z scale
      eps = 1.d-7
      mean = 0.d0
      qs = dgauin(1-tr)
      do 300 i=1,k

c       probability of crossing upper boundary by pt(i)
c
         pn(i) = alpha1(pt(i),alpha)
         pd = pn(i) - pn(i-1)
         var = pt(i) - pt(i-1)
         flag2=0

         if(pd .le. 0.d0) then
            yn(i) = 5.d0
c          nyn(i) = -5.d0
            goto 175
         endif

         if(i .eq. 1) then
            yn(i) = dgauin(1.d0-pd)*dsqrt(var) + mean
            goto 175
         endif

c       integrate fcab from infinity to yn(i)

         sum = 0.d0
         dnp = 5.d0
         fold = fcab(dnp,i,mean,var)
 100     dnp = dnp - h
         fnew = fcab(dnp,i,mean,var)
         temp = h * (fold+fnew) / 2.d0
         sum = sum + temp

c       check if upper limit is large enough

         if(sum .ge. pd) goto 150
         if (dnp.lt.0) then
            ierr=1
            return
         endif
         fold = fnew
         goto 100

c       initialize for falsi3

 150     dnp = dnp + h
         sum = sum - temp
         a = 0.d0
         b = h
         flag = -1

c call falsi3 for interpolation

         yn(i) = falsi3(flag,dnp,fold,i,a,b,sum,temp,pd,mean,var)

c check if greater than truncated boundary
 175     temp=qs*sqrt(pt(i))
         if (yn(i).gt.temp) then
            yn(i)=temp
            flag2=1
         endif

c determine the number of grids from -infinity to yn(i)

c        nyn(i) = -5.d0
         ln(i) = (yn(i) - nyn(i)) / h

         if(i .eq. k) goto 399

c       compute fcab for later use

         lnp2 = ln(i) + 2
         if (lnp2.gt.500) then
            ierr=1
            return
         endif
         dnp = nyn(i)
         do 200 j=1,lnp2
            if(j .eq. lnp2) dnp = yn(i)
            fn(j,i) = fcab(dnp,i,mean,var)
            dnp = dnp + h
 200     continue
c compute true rejection probability
c fn(j,i) is the joint density of the stopping time (i) and the value of the 
c statistic.  Rejects at i if stat>=yn(i), so Prob(stopping time = i)
c = \int_{yn(i)}^{\infty} f(j,i) dj 
 399     if (i.eq.1.and.flag2.eq.1) then
            pn(i)=tr
         else if (yn(i).ge.5) then
            pn(i)=0
         else
            sum=0
            dnp=yn(i)
            fold = fcab(dnp,i,mean,var)
            lnp2=(5-dnp)/h
            do 288 j=1,lnp2
               dnp=dnp+h
               fnew = fcab(dnp,i,mean,var)
               temp = h * (fold+fnew) / 2.d0
               sum = sum + temp
               fold=fnew
 288        continue
            pn(i)=sum+pn(i-1)
c	temp=(fn(1,i)+fn(ln(i)+1,i))/2
c	do 208 j=2,ln(i)
c	   temp=temp+fn(j,i)
c 208	continue
c	temp=temp*h+(fn(ln(i)+1,i)+fn(ln(i)+2,i))*
c     1    (yn(i)-nyn(i)-h*ln(i))/2
c	dnp=alpha1(pt(i),alpha)
c        else
c           pn(i)=1-temp
         endif
 300  continue
 400  return
      end

      real*8 function falsi3(flag,dnp,fold,k,a,b,sum,temp,pd,mean,var)

c linear interpolation procedure of the modified regula falsi method

c   calling sequence

      integer k,flag
      real*8 dnp,pd,a,b,sum,temp,fold,mean,var

c   local definitions

      real*8 fa,fb,dnpw,w,fw,fnew,signfa,prevfw,eps,zero
      real*8 fcab

      eps = 1.d-7
      zero = 0.d0
      fa = sum - pd
      fb = fa + temp
      signfa = dsign(1.d0,fa)
      
c a and b are such that fa*fb < 0.

      w = a
      fw = fa
 1    w = (fa*b - fb*a) / (fa - fb)
      prevfw = dsign(1.d0,fw)
      dnpw = dnp + flag * w
      fnew = fcab(dnpw,k,mean,var)
      temp = w * (fold+fnew) / 2.d0
      fw = sum + temp - pd
      if(dabs(fw) .lt. eps) goto 2
      if(signfa*fw .lt. zero) then
         b = w
         fb = fw
         if(fw*prevfw .gt. zero) fa = fa / 2.d0
      else
         a = w
         fa = fw
         if(fw*prevfw .gt. zero) fb = fb / 2.d0
      endif
      goto 1

 2    falsi3 = dnpw

      return
      end
