      double precision function schwarz_shell(ish, jsh)
C$Id: schwarz_sh.F 19707 2010-10-29 17:59:36Z d3y133 $
      implicit none
#include "cschwarz.fh"
#include "mafdecls.fh"
      integer ish, jsh
      double precision value
#include "itri.fh"
c
      call comp4_extract(int_mb(k_csh), itri(ish,jsh), value, nb_per_i)
      schwarz_shell = value
c
      end
