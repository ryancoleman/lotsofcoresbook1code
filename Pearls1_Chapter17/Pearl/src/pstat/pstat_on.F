      Subroutine PStat_On( Handle )
C$Id: pstat_on.F 19708 2010-10-29 18:04:21Z d3y133 $
      Implicit NONE
#include "errquit.fh"
      Integer Handle
C
#include "pstatP.fh"
#include "util.fh"
C
      double precision perfm_flop
      external perfm_flop
      Double precision Zero
      Parameter (Zero = 0.0d0)
      logical status
C
      status = handle.gt.0 .and. handle.le.pstat_max_items
      if (status) status = pstat_active(handle)
CHF   if (.not. status) Call ErrQuit
CHF  $     ('pstat_on: Timer not assigned', Handle)
      if (.not. status) return
c
      If ( PStat_Unmatched_Starts( Handle ) .ne. 0  ) then
         Call ErrQuit('pstat_on: Timer is already on', Handle,
     &       UNKNOWN_ERR)
      Else
         PStat_Unmatched_Starts( Handle ) = 1
      EndIf
C
      If ( PStat_Functions( Handle, PStat_QStat) ) then
         PStat_Start_Times( Handle, 1 ) = util_cpusec()
         PStat_Start_Times( Handle, 2 ) = util_wallsec()
         pstat_start_flop(handle) = perfm_flop()
      EndIf
C
      Return
      End
