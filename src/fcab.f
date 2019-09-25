      real*8 function fcab (x,k,mean,varian)

c this function computes recursively defined density functions of (S,N).

c calling sequence
c      x := value at which density function is evaluated
c      k := number of analysis
c   mean := mean of the Gaussian process
c varian := variance of the Gaussian process

      integer k
      real*8 x,mean,varian

c local definitions

      integer j,jln,km1
      real*8 a,temp,stdv
      real*8 phi

c common definitions

      integer ln(30)
      real*8 fn(500,30),nyn(30),yn(30),h
      common /amr/ fn,nyn,yn,h,ln

c common definition with ones and fcab
c    h   := grid size for numerical integration by trapezoidal rule
c    ln  := number of grids
c    fn  := coordinate of the recursively defined joint density
c    yn  := upper limit of integration (boundary)
c    nyn := lower limit of integration (boundary)

      stdv = dsqrt(varian)
      if(k .eq. 1) then
         temp = phi((x-mean)/stdv)/stdv
      else
         km1 = k - 1
         jln = ln(km1)
         a = nyn(km1)
         temp = fn(1,km1) * phi((x-mean-a)/stdv)/stdv

c integrate fcab(*,km1,mean,varian) from nyn(km1) to nyn(km1) + jln*h

         do 100 j=2,jln
            a = a + h
            temp = temp + 2.d0*fn(j,km1)*phi((x-mean-a)/stdv)/stdv
 100     continue

c area at the end of the grid

         a = a + h
         temp = temp + fn(jln+1,km1) * phi((x-mean-a)/stdv)/stdv
         temp = h / 2.d0 * temp

c area of the last interval of length less than h

         temp = temp + (yn(km1)-a) / 2.d0
     +        * (fn(jln+1,km1)*phi((x-mean-a)/stdv)/stdv
     +        +fn(jln+2,km1)*phi((x-mean-yn(km1))/stdv)/stdv)
      endif

      fcab = temp
      return
      end
