      real*8 function phi(x)
        
c This function subprogram computes the probability density function of
c the standard normal distribution.

      real*8 x

      real*8 pi,sq2pi

      parameter (pi=3.1415926535898d0)

      sq2pi=dsqrt(2.d0*pi)

      phi=dexp(-x*x/2.d0)/sq2pi

      return
      end

      real*8 function dgaudf(z)

c This function subprogram computes the distribution function of the
c standrad normal distribution.  The function subprogram `phi' is
c used to evaluate the standard noraml density function.

c This routine is more accurate than the IMSL routine `dnordf'.

c The original version had been written by Russ Spry at Biostatistics
c Center, University of Wisconsin, Madison.  It was last modified on
c November 24, 1989, by Kyungmann Kim at Division of Biostatistics,
c Dana Farber Cancer Center, Boston, Massachusetts, USA.

      real*8 z

      integer k,n
      real*8 a,c,s,u
      real*8 phi

      parameter (n=16,a=2.8d0)

      if (dabs(z).lt.a) then
         s=1.d0
         u=1.d0
         do 1 k=1,n
            u=u*z*z/(2*k+1)
 1          s=s+u
         dgaudf=5.d-1+z*s*phi(z)
      else
         c=0.d0
         do 2 k=n,1,-1
 2       c=k/(z+c)
         c=1.d0/(z+c)*phi(z)
         if (z.gt.0.d0) then
            dgaudf=1.d0-c
         else
            dgaudf=-c
         endif
      endif

      return
      end

      real*8 function dgauin (p)

c This function subprogram computes the quantile of the standard normal
c distribution.  The function subprogram `phi' is used to evaluate
c the standard noraml density function.

c This routine is more accurate than the IMSL routine `dnorin'.

c The original version had been written by Russ Spry at Biostatistics
c Center, University of Wisconsin, Madison.  It was last modified on
c November 24, 1989, by Kyungmann Kim at Division of Biostatistics,
c Dana Farber Cancer Center, Boston, Massachusetts, USA.

      real*8 p

      real*8 eps,errrel,x0,x1
      real*8 dgaudf,phi

      data eps /1.d-12/

c use the Newton-Raphson method to find the solution to p=dgaudf(z)

      if (p .eq. 0.5d0) then
         dgauin=0.d0
      else
         x0=0.d0
 1       x1=x0-(dgaudf(x0)-p)/phi(x0)
         errrel=(x0-x1)/x1
         if (dabs(errrel).gt.eps) then
            x0=x1
            goto 1
         endif
         dgauin=x1
      endif

      return
      end


      real*8 function alpha1 (t,pd)

c This use function corresponds to the alpha 1 star in Lan & DeMets(1983).

      real*8 t,pd

c anordf and anorin are function subroutines needed for this alpha star.

      real*8 dgaudf,dgauin

      alpha1 = 2.d0 * (1.d0 - dgaudf(dgauin(1.d0-pd/2.d0)/dsqrt(t)))
      return
      end


      real*8 function alpha2 (t,pd)
c
c This use function corresponds to alpha 2 star in Lan & DeMets (1983).
c
      real*8 t,pd,e
c
      e = 2.718281828459d0
      alpha2 = pd * dlog(1.d0 + (e-1.d0)*t)
      return
      end


      real*8 function alpha3 (t,pd)
c
c This use function corresponds to alpha 3 star in Lan & DeMets (1983).
c
      real*8 t,pd
c
      alpha3 = pd * t
      return
      end


      real*8 function alpha4 (t,pd)

c This use function corresponds to the alpha 4 star in Kim & DeMets (1987).

      real*8 t,pd
      
      alpha4 = pd * (t ** 1.5)

      return
      end


      real*8 function alpha5 (t,pd)

c This use function corresponds to the alpha 5 star in Kim & DeMets (1987).

      real*8 t,pd

      alpha5 = pd * t * t

      return
      end

      real*8 function nfa(rct,ll,A,jj,pi)

c calculate the number of failures during the accrual period

      integer j,jj
      real*8 rct,ll(30),A,pi(30)

      nfa = 0.d0
      do 1 j=1,jj
         nfa = nfa + A*pi(j)*(rct-(1.d0-dexp(-ll(j)*rct))/ll(j))
 1    continue
      return
      end

      real*8 function nff(rct,ll,A,sa,jj,pi)

c calculate the number of failures during the ffolow-up period

      integer j,jj
      real*8 rct,ll(30),A,sa,pi(30)

      nff = 0.d0
      do 1 j=1,jj
         nff = nff 
     +        + A*pi(j)*(sa-(dexp(ll(j)*sa)-1.d0)/dexp(ll(j)*rct)/ll(j))
 1    continue
      return
      end
