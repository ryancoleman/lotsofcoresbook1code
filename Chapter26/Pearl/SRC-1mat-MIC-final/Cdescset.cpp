#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cdescset(int *desc, int m, int n,int mb,int nb,
               int irsrc,int icsrc,int ictxt, int lld )
{
  /*
   * Like Cdescinit but no checking
   */

  desc[DTYPE_] = BLOCK_CYCLIC_2D;
  desc[M_]     = m;
  desc[N_]     = n;
  desc[MB_]    = mb;
  desc[NB_]    = nb;
  desc[RSRC_]  = irsrc;
  desc[CSRC_]  = icsrc;
  desc[CTXT_]  = ictxt;
  desc[LLD_]   = lld;

  return;
}

