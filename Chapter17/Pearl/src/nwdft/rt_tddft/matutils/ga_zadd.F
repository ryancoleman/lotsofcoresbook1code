C
C     Hack for GA 4-3 which lacks ga_zadd(). GNUmakefile will only
C     include if OLD_GA.
C
      subroutine ga_zadd(za, g_a, zb, g_b, g_c)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer g_a, g_b, g_c
      double complex za, zb
      
      call ga_add (za, g_a, zb, g_b, g_c)

      end subroutine
c $Id: ga_zadd.F 24344 2013-06-22 06:15:21Z d3y133 $
