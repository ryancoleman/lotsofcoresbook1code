#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cdescinit(int *desc, int m, int n,int mb,int nb,
               int irsrc,int icsrc,int ictxt, int lld, int *info)
{
  int isok;
  int Locp;
  int nprow,npcol,   myprow,mypcol;

  const int idebug  = 1;

  *info = 0;

  Cblacs_gridinfo( ictxt, &nprow, &npcol, &myprow, &mypcol );

  isok = (nprow >= 1) && (npcol >= 1);;
  if (!isok) {
    *info = -8;
    if (idebug >= 1) {
      printf("Cdescinit: nprow %d npcol %d \n",
                         nprow,   npcol );
    };
    return;
  };

  isok = (0 <= myprow) && (myprow < nprow);
  if (!isok) {
    *info = -8;
    if (idebug >= 1) {
      printf("Cdescinit: myprow %d nprow %d \n",
                         myprow,   nprow );
    };
    return;
  };

  isok = (0 <= mypcol) && (mypcol < npcol);
  if (!isok) {
    *info = -8;
    if (idebug >= 1) {
      printf("Cdescinit: mypcol %d npcol %d \n",
                         mypcol,   npcol);
    };
    return;
  };

  if (m <= 0) {
    *info = -2;
    return;
  };

  if (n <= 0) {
    *info = -3;
    return;
  };
  if (mb < 1) {
    *info = -4;
    return;
  };
  if (nb < 1) {
    *info = -5;
    return;
  };

  isok = (0 <= irsrc) && (irsrc < nprow);
  if (!isok) {
    *info = -6;
    return;
  };

  isok = (0 <= icsrc) && (icsrc < npcol);
  if (!isok) {
    *info = -7;
    return;
  };


  Locp = Cnumroc(m,mb,myprow,irsrc,nprow);
  isok = (lld >= MAX(1,Locp));
  if (!isok) {
    *info =  -9;
    if (idebug >= 1) {
       printf("Cdescinit: m %d mb %d myprow %d irsrc %d nprow %d \n",
                          m,   mb,   myprow,   irsrc,   nprow );
       printf("Cdescinit: Locp %d lld %d \n",
                          Locp,   lld );
       };
    return;
    };


    desc[ DTYPE_ ] = BLOCK_CYCLIC_2D;
    desc[ M_     ] = m;
    desc[ N_     ] = n;
    desc[ MB_    ] = mb;
    desc[ NB_    ] = nb;
    desc[ RSRC_  ] = irsrc;
    desc[ CSRC_  ] = icsrc;
    desc[ CTXT_  ] = ictxt;
    desc[ LLD_   ] = MAX( lld, MAX(1, Locp) );

    return;
}
