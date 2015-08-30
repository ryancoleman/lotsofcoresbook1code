      subroutine scatter(n,a,indx,b)
*
* $Id: scatter.f 19707 2010-10-29 17:59:36Z d3y133 $
*
      integer n, indx(n)
      double precision a(*), b(n)
      integer i
      
      do i=1,n
        a(indx(i)) = b(i)
      enddo
      return
      end
