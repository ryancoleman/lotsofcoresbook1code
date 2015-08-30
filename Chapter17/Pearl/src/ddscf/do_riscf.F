      subroutine do_riscf(doit)
C$Id: do_riscf.F 19707 2010-10-29 17:59:36Z d3y133 $
      implicit none
#include "cfock.fh"
      logical doit
      
      rifock = doit

      return
      end
