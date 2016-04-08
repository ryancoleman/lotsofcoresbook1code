      Subroutine dft_scatter(n,a,index,b)

C$Id: dft_scatter.F 19696 2010-10-29 16:53:42Z d3y133 $
      implicit none
      
      integer n
      integer index(n)
      double precision a(*),b(n)
      integer i

c*******************************************************************************

      do 10 i = 1,n
       a(index(i)) = b(i)
   10 continue

      end
