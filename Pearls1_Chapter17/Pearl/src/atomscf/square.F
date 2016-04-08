      subroutine square(r,a,mrowr,n)
C$Id: square.F 19708 2010-10-29 18:04:21Z d3y133 $
      implicit none
      integer mrowr, n
      real *8  r(*),a(*)
c
c... convert triangle a to square r
c
      integer i,j
      integer iii,jj,jjj,k
      iii=0
      jjj=1
      k=1
      do  i=1,n
         jj=jjj
         do   j=1,i
           r(iii+j)=a(k)
           r(jj)=a(k)
           k=k+1
           jj=jj+mrowr
         enddo
         iii=iii+mrowr
         jjj=jjj+1
      enddo
      return
      end
