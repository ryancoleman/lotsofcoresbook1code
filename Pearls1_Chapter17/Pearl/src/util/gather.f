
      subroutine gather(n,a,b,indx)
*
* $Id: gather.f 19707 2010-10-29 17:59:36Z d3y133 $
*
      integer n, indx(n)
      double precision a(n),b(*)
      integer i
      
      do i=1,n
        a(i) = b(indx(i))
      enddo
      return
      end

      






