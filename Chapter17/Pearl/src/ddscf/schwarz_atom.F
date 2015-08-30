      double precision function schwarz_atom(iat, jat)
C$Id: schwarz_atom.F 19707 2010-10-29 17:59:36Z d3y133 $
      implicit none
#include "cschwarz.fh"
#include "mafdecls.fh"
      integer iat, jat
      double precision value
#include "itri.fh"
c
      call comp4_extract(int_mb(k_cat), itri(iat,jat), value, nb_per_i)
      schwarz_atom = value
c
      end
