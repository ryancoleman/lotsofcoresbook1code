c $Id: hnd_terminate.F 19696 2010-10-29 16:53:42Z d3y133 $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : termination routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function hnd_terminate()
      implicit none
#include "apiP.fh"
c
c     at this point nothing to do but just in case for later
c
      hnd_terminate = .true.
      return
      end
