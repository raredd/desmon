      subroutine squse(alpha,alphal,kk,pt,use,usel,tr,trl,eta,uz,
     $     uzl,ierr)
c alpha=one sided sig level
c use=type use function 1=OF, 2=Pocock, 3=linear, 4=1^1.5, 5=t^2
c kk # analysis times; pt information times (pt(0)=0 is not an analsis time)
c eta=brownian motion drift paramter for alternative (needed for rci)
      double precision alpha,pt(0:kk),alphal,tr,trl,eta
      integer kk,use,usel,ierr

c  ********************************************************************
c    this routine calculates the group sequential boundaries, nominal
c               significance, and rejection probabilities
c  ********************************************************************

      integer k
      real*8 pn(0:kk),uz(kk),temp,uzl(kk)

c common definition with ones and fcab
c    h   := grid size for numerical integration by trapezoidal rule
c    ln  := number of grids
c    fn  := coordinate of the recursively defined joint density
c    yn  := upper limit of integration (boundary)
c    nyn := lower limit of integration (boundary)

      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

      h = 0.05d0
      pt(0) = 0.d0
      pn(0) = 0.d0

c call ones to generate one-sided boundaries - first set lower limit
c for one-sided boundary
c generate one-sided group sequential boundary using use function for rci

         do 10 k=1,kk
            nyn(k)=-5
 10      continue 
      if (alphal.gt.0) then
         call ones(kk,alphal,pt,pn,usel,trl,ierr)
         if (ierr.gt.0) return
c compute approx rci lower boundary
         do 810 k=1,kk
            nyn(k)=-yn(k)+eta*pt(k)
 810     continue 
      endif

c compute upper boundary

      call ones(kk,alpha,pt,pn,use,tr,ierr)
      if (ierr.gt.0) return

c standardize boundaries

      do 30 k=1,kk
         temp = dsqrt(pt(k))
         uz(k) = yn(k)/temp
         uzl(k) = nyn(k)/temp
 30   continue

c$$$
c$$$      do 10 k=1,kk
c$$$         nyn(k)=-5
c$$$ 10   continue 
c$$$
c$$$      call ones(kk,alpha,pt,pn,use,temp,tr)
c$$$
c$$$c determine the boundaries in standardized form
c$$$
c$$$      do 40 k=1,kk
c$$$        uz(k) =  yn(k) / dsqrt(pt(k))
c$$$40    continue
c$$$
c$$$c k is the analysis number, pt(k) the information time uz(k) the boundary
c$$$c dgaudf(-uz(k)) the nominal sig level, pn(k) is the significance level
c$$$c used at pt(k) for this use function
c$$$c      write(iounit,998) (k,pt(k),uz(k),dgaudf(-uz(k)),
c$$$c     +pn(k)-pn(k-1),k=1,kk)
      return
      end

