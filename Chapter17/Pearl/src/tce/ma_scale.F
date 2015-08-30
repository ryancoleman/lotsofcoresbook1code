      subroutine ma_scale(a,n,scalar)
c
c $Id: ma_scale.F 19706 2010-10-29 17:52:31Z d3y133 $
c
      implicit none
#include "mafdecls.fh"
#include "tce.fh"
      integer i,n
      double precision a(n)
      double precision scalar
c
      do i = 1,n
        a(i) = a(i) * scalar     
      enddo
      return
      end
