      Logical Function PStat_Terminate()
C$Id: pstat_term.F 19708 2010-10-29 18:04:21Z d3y133 $
      Implicit NONE
C
#include "pstatP.fh"
C
C     Nothing to do here, really, but let us "unassign" all of the items
C
      integer i
c
      do i = 1, pstat_max_items
         pstat_active(i) = .false.
      enddo
C
      PStat_Terminate = .TRUE.
      Return
      End
