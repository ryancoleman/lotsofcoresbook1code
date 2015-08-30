c
      LOGICAL FUNCTION ACOLO_1P(p1b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER p1b
      INTEGER iocc_v
      acolo_1p = .false.
      iocc_v=int_mb(k_active+p1b-1)
      if(iocc_v.ge.numact-2) acolo_1p=.true.
      RETURN
      END
c $Id: acolo_1p.F 21176 2011-10-10 06:35:49Z d3y133 $
