#include "ooclu.h"



#ifdef __cplusplus
extern "C" 
#endif
void setup_desc(int m,int n,int ia,int ja,int *descA,      
                int *isize, int *descB )
{
/*
 * estimate the amount of storage and setup descB
 */

 const int use_round_up_lld = FALSE;
 const int use_actual_lld = FALSE;
 const int align_size = 16;
 const int use_round_up_mn = FALSE;

 int  mr, nr;
 int  ictxt = descA[CTXT_];
 int  lld, Locp, Locq, rsrc, csrc;
 int  Locp_actual, Locq_actual;
 int  nprow, npcol, myprow, mypcol;
 int  mb,nb,   irsrc,icsrc;

 Cblacs_gridinfo( ictxt, &nprow, &npcol, &myprow, &mypcol );


 /* 
  * make B(1,1) aligned to A(ia,ja)
  */

 mb = descA[MB_];
 nb = descA[NB_];

 /*
  * find out processor (rsrc,csrc) that owns A(ia,ja)
  */
 irsrc = descA[RSRC_];
 icsrc = descA[CSRC_];

 /*
  * check for replicated storage format
  */
 if (irsrc == -1) {
   rsrc = -1;
   }
 else {
   rsrc = Cindxg2p( ia, mb, myprow, irsrc, nprow );
 };

 if (icsrc == -1) {
   csrc = -1;
   }
 else {
 csrc = Cindxg2p( ja, nb, mypcol, icsrc, npcol );
 };



 /*
  * Over estimate amount of local storage
  * use same amount of storage on all processors
  */

 mr = m;
 nr = n;
 if (use_round_up_mn) {
   /*
    * round up to next block size
    */
   mr = ((m + (mb-1))/mb)*mb;
   nr = ((n + (nb-1))/nb)*nb;
 };


 Locp = Cnumroc( mr, mb, 0, 0, nprow );
 Locp_actual = Cnumroc( mr, mb, myprow, rsrc, nprow);

 Locq = Cnumroc( nr, nb, 0, 0, npcol );
 Locq_actual = Cnumroc( nr,nb, mypcol, csrc,npcol);


 lld = MAX(1,Locp);
 if (use_actual_lld) {
   lld = MAX(1, Locp_actual);
 };

 if (use_round_up_lld) {
   lld = ((lld + (align_size-1))/align_size)*align_size;
 };

 *isize = MAX(1,lld) * MAX(1,Locq);

 descB[DTYPE_]  = descA[DTYPE_];
 descB[M_]      = mr;
 descB[N_]      = nr;
 descB[MB_]     = descA[MB_];
 descB[NB_]     = descA[NB_];
 descB[RSRC_]   = rsrc;
 descB[CSRC_]   = csrc;
 descB[CTXT_]   = descA[CTXT_];
 descB[LLD_]    = lld;

 return;
}

