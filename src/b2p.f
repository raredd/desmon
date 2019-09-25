      subroutine pfishr(alpha,p1,p2,n1,n2,d,cd,d1,d2,cd1,cd2,cum,cumr)
      double precision d(0:2000),cd(0:2000),nd1,md,alpha
      double precision d1(0:1000),d2(0:1000),cd1(0:1000),cd2(0:1000)
      double precision nd2,p1,p2,al2,cumr,t1,cum
      integer n1,n2,n,i,j,k,k1,k2
      al2=1-alpha
      nd1=n1
      nd2=n2
      call bin(nd1,p1,d1,cd1)
      call bin(nd2,p2,d2,cd2)
      n=n1+n2-1
      cumr=d1(0)*d2(0)*alpha
      cumr=cumr+d1(n1)*d2(n2)*alpha
      cum=0.d0
      do 10 i=1,n
         md=dfloat(i)
         call hyp(nd1,nd2,md,d,cd)
         k1=min0(i,n1)
         k=k1
         if (cd(k1-1).lt.al2) go to 11
         t1=d1(k1)*d2(i-k1)
         cum=cum+t1
         cumr=cumr+t1
         k1=k1-1
         k=k1
         k2=max0(0,i-n2)
         do 12 j=k1,k2,-1
            if (cd(k-1).lt.al2) go to 11
            t1=d1(k)*d2(i-k)
            cum=cum+t1
            cumr=cumr+t1
            k=k-1
 12      continue
 11      if (d(k).eq.0.d0) go to 10
         cumr=cumr+(cd(k)-al2)*d1(k)*d2(i-k)/d(k)
 10   continue
 999  return
      end

      subroutine bin(n1,p1,d,cd)
      double precision d(0:1000),cd(0:1000),p1,q1,n1,o1,cum
      integer i,n
      q1=1-p1
      o1=dlog(p1/q1)
      n=idint(n1+.5)
      do 10 i=0,n
      d(i)=0.d0
  10  cd(i)=0.d0
      cum=n*dlog(q1)
      if (cum.lt.-200) go to 8
      d(0)=dexp(cum)
      go to 9
   8  d(0)=0.d0
   9  cd(0)=d(0)
      do 15 i=1,n
      cum=cum+o1+dlog(n-dfloat(i)+1.d0)-dlog(dfloat(i))
      if (cum.lt.-200.d0) go to 26
      d(i)=dexp(cum)
      cd(i)=cd(i-1)+d(i)
      go to 15
  26  d(i)=0.d0
      cd(i)=cd(i-1)
  15  continue
      return
      end

      subroutine hyp(n1d,n2d,md,d,cd)
      double precision d(0:2000),cd(0:2000),cum,n1d,n2d,md,id,id1
      integer n1,n2,m,ll,lu,i,n,j,l,l2
      n1=idint(n1d+.5)
      n2=idint(n2d+.5)
      m=idint(md+.5)
      ll=max0(0,m-n2)
      lu=min0(n1,m)
      n=n1+n2
      do 10 i=0,n
      d(i)=0
  10  cd(i)=0
      cum=0
      do 12 i=2,n
  12  cum=cum-dlog(dfloat(i))
      if (m.lt.2) go to 31
      do 13 i=2,m
  13  cum=cum+dlog(dfloat(i))
  31  j=n-m
      if (j.lt.2) go to 32
      do 14 i=2,j
  14  cum=cum+dlog(dfloat(i))
  32  l=n1
      if (ll.eq.0) l=n2
      l2=ll
      if (ll.eq.0) l2=m
      do 15 i=2,l
  15  cum=cum+dlog(dfloat(i))
      if (l2.lt.2) go to 33
      do 16 i=2,l2
  16  cum=cum-dlog(dfloat(i))
  33  j=l-l2
      if (j.lt.2) go to 34
      do 17 i=2,j
  17  cum=cum-dlog(dfloat(i))
  34  if (cum.lt.-200.d0) go to 53
      d(ll)=dexp(cum)
      cd(ll)=d(ll)
      go to 54
  53  d(ll)=0.d0
      cd(ll)=0.d0
  54  ll=ll+1
      do 20 i=ll,lu
      id=dfloat(i-1)
      id1=dfloat(i)
      cum=cum+dlog((n1-id)*(md-id)/(id1*(n2-md+id1)))
      if (cum.lt.-200.d0) go to 26
      d(i)=dexp(cum)
      cd(i)=cd(i-1)+d(i)
      go to 20
  26  d(i)=0.d0
      cd(i)=cd(i-1)
  20  continue
      return
      end
