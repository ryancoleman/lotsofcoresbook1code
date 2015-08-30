      subroutine dft_get2eri(CD_bas_han, g_2ceri,oskel)
      implicit none
      integer CD_bas_han, g_2ceri(*)
      logical oskel
      integer maxg, mscratch_2e2c
      external int_2e2c
c
      call int_mem_2e2c(maxg, mscratch_2e2c)
      call ga_get2eri(cd_bas_han, g_2ceri,oskel,maxg, mscratch_2e2c, 
     &   int_2e2c)
c
      return
      end
c $Id: dft_get2eri.F 21176 2011-10-10 06:35:49Z d3y133 $
