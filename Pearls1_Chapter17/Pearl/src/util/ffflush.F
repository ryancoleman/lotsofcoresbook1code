      subroutine ffflush(unit)
C$Id: ffflush.F 19707 2010-10-29 17:59:36Z d3y133 $
      integer unit              ! [input]
c
c     Flush outputto this fortran unit
c
      call util_flush(unit)
c
      end
