
#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cinfog1l( int gindx, int nb, int nprocs, int myroc, int isrcproc, 
               int *lindx_out, int *rocsrc_out )
{
  int gcpy, iblk, lindx, rocsrc;

  gcpy = gindx-1;
  iblk = gcpy / nb;
  rocsrc = MOD( iblk + isrcproc, nprocs );
  lindx = ( iblk / nprocs + 1 ) * nb + 1;

  if( MOD(myroc+nprocs-isrcproc,nprocs) >= MOD(iblk, nprocs) ) {
     if( myroc == rocsrc ) {
       lindx = lindx + MOD( gcpy, nb );
       };
      lindx = lindx - nb;
  };
 
  *lindx_out = lindx;
  *rocsrc_out = rocsrc;

  return;
}
