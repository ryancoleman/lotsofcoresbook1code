      Subroutine dft_gather(n,a,b,index,idiff)

C$Id: dft_gather.F 19696 2010-10-29 16:53:42Z d3y133 $
      implicit none
      
      integer n
      double precision a(n),b(*)
      integer index(n),idiff
      integer i

c****************************************************************************

      do 10 i = 1,n
       a(i) = b(index(i)-idiff)
   10 continue

      end
