c 6 or 7 works well here
#define UNROLLING 6
c
      subroutine ma_transpose(a,n)
c
c $Id: ma_transpose.F 19706 2010-10-29 17:52:31Z d3y133 $
c
      implicit none
c#include "mafdecls.fh"
c#include "tce.fh"
      integer i,j,n
      double precision a(n,n)
      double precision swap
c
!DEC$ prefetch a
!DEC$ unroll(UNROLLING)
!DEC$ vector always
      do i = 1,n
!DEC$ prefetch a
!DEC$ unroll(UNROLLING)
!DEC$ vector always
        do j = i,n
          swap = a(i,j)
          a(i,j) = a(j,i)      
          a(j,i) = swap
        enddo
      enddo
      return
      end
