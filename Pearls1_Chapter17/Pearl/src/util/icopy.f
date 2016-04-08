      Subroutine icopy(n,src,ix,dest,iy)
*
* $Id: icopy.f 19707 2010-10-29 17:59:36Z d3y133 $
*
      implicit none
      integer n, src(*), ix, dest(*), iy
      integer i,ii,jj

      ii = 1
      jj = 1
      do i=1,n
        dest(jj) = src(ii)
        ii = ii + ix
        jj = jj + iy
      enddo
      return
      end

