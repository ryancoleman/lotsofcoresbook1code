      subroutine riscf_dovecs ( doit, g_vec )
C$Id: riscf_dovecs.F 19707 2010-10-29 17:59:36Z d3y133 $
C     switches on use of MO vectors (Ken Dyalls proposal) for RI-SCF
      implicit none
#include "cfock.fh"
      logical doit
      integer g_vec

      dovecs = doit
      g_vecs = g_vec

      return
      end
