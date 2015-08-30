      double precision function dabssum(n, a)
c $Id: dabssum.F 19707 2010-10-29 17:59:36Z d3y133 $
      implicit none
c
      integer n, i
      double precision a(n)
c
      dabssum = 0.0d0
      do i = 1, n
         dabssum = dabssum + abs(a(i))
      enddo
c
      end
