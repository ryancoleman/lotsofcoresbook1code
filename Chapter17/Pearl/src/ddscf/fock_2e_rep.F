      subroutine fock_2e_report(kat, natoms)
*
* $Id: fock_2e_rep.F 19707 2010-10-29 17:59:36Z d3y133 $
*
      implicit none
#include "cfock.fh"
#include "global.fh"
c
c     Print out a reassuring message every few iterations
c     of the k atom loop during fock matrix construction.
c
      integer kat, natoms
      character*26 date
      logical reportit
c
      reportit = ((mod(kat,10).eq.0) .or. (kat.eq.1) .or.
     $     (kat.eq.natoms)) .and. (ga_nodeid().eq.0)
      if (reportit) then
         call util_date(date)
         write(6,1) kat, date
 1       format(' fock_2e : katom=',i4,', date=',a)
         call util_flush(6)
      endif
c
      end

