#include "ooclu.h"












#ifdef __cplusplus
extern "C" 
#endif
void local_extent_info(int m,int n,int ia,int ja,int *descA, 
                  int *msize, int *nsize, 
                  int *lr1, int *lc1, int *lr2, int *lc2 )
{
/*
 * compute the local extent of the distributed submatrix
 * (lr1,lc1)  is the begining local Fortran array index
 * (lr2,lc2)  is the ending local Fortran array index
 * The LOCAL  part of the distributed matrix on this processor is
 * A( lr1:lr2, lc1:lc2 )
 *
 *
 * Conceptually number of local rows is  msize = (lr2-lr1+1)
 * and  number of local columns is nsize = (lc2-lc1+1)
 *
 * However, if msize <= 0, or nsize <= 0, then (lr1,lc1), (lr2,lc2)
 * may not be valid
 */
  int icontxt, nprow,npcol,  myprow,mypcol;
  int mb,nb,rsrc,csrc;
  int jlr1,jlr2,  jlc1,jlc2;
  int ia_first,ia_last, ja_first, ja_last, iproc;

  const int use_infog2l = FALSE;

  mb = descA[MB_];
  nb = descA[NB_];
  rsrc = descA[RSRC_];
  csrc = descA[CSRC_];
  icontxt = descA[CTXT_];

  Cblacs_gridinfo( icontxt, &nprow,&npcol,   &myprow,&mypcol);



  jlr1 = -1;
  jlr2 = -1;

  ia_first = Ciafirst( ia, mb,myprow,rsrc,nprow);
  ia_last = Cialast( ia+m-1,mb,myprow,rsrc,nprow);

  if ( (ia <= ia_first) && (ia_first <= ia + m-1)) {
    Cinfog1l( ia_first, mb, nprow, myprow, rsrc,
           &jlr1,  &iproc );
    assert( iproc == myprow );
  };

  if ((ia <= ia_last) && (ia_last <= ia+m-1)) {
    Cinfog1l( ia_last, mb, nprow, myprow, rsrc,
          &jlr2, &iproc );
    assert( iproc == myprow);
  };



  jlc1 = -1;
  jlc2 = -1;

  ja_first = Ciafirst( ja,nb,mypcol,csrc,npcol);
  ja_last = Cialast( ja+n-1,nb,mypcol,csrc,npcol);


  if ((ja <= ja_first) && (ja_first <= ja + n-1)) {
    Cinfog1l( ja_first, nb, npcol, mypcol, csrc,
           &jlc1, &iproc );
    assert( iproc == mypcol );
  };

  if ((ja <= ja_first) && (ja_first <= ja + n-1)) {
    Cinfog1l( ja_last, nb, npcol, mypcol, csrc,
           &jlc2, &iproc );
    assert( iproc == mypcol );
  };


  *msize = 0;
  if ((1 <= jlr1) && (jlr1 <= jlr2)) {
    *msize = jlr2 - jlr1 + 1;
  };
  *lr1 = jlr1;
  *lr2 = jlr2;

  *nsize = 0;
  if ((1 <= jlc1) && (jlc1 <= jlc2)) {
    *nsize = jlc2 - jlc1 + 1;
  };
  *lc1 = jlc1;
  *lc2 = jlc2;

  return;
}




#ifdef __cplusplus
extern "C" 
#endif
void local_extent_loop(int m,int n,int ia,int ja,int *descA, 
                  int *msize, int *nsize, 
                  int *lr1, int *lc1, int *lr2, int *lc2 )
{
/*
 * compute the local extent of the distributed submatrix
 * (lr1,lc1)  is the begining local Fortran array index
 * (lr2,lc2)  is the ending local Fortran array index
 * The LOCAL  part of the distributed matrix on this processor is
 * A( lr1:lr2, lc1:lc2 )
 *
 *
 * Conceptually number of local rows is  msize = (lr2-lr1+1)
 * and  number of local columns is nsize = (lc2-lc1+1)
 *
 * However, if msize <= 0, or nsize <= 0, then (lr1,lc1), (lr2,lc2)
 * may not be valid
 */

 const int idebug = 0;

 int i,j,iproc;

 int nprow,npcol,myprow,mypcol;
 int lrindx,lcindx;
 int lrindx1,lcindx1, lrindx2,lcindx2, rsrc,csrc, rproc,cproc;
 int mb,nb;
 int iia,jja;

 int is_replicated;
 int isok;
 int is_found;

 Cblacs_gridinfo( descA[CTXT_], &nprow,&npcol,  &myprow,&mypcol);
 assert( nprow >= 1);
 assert( npcol >= 1);
 assert( (0 <= myprow) && (myprow < nprow));
 assert( (0 <= mypcol) && (mypcol < npcol));

 rsrc = descA[RSRC_];
 csrc = descA[CSRC_];

 if (rsrc == -1) { rsrc = myprow; };
 if (csrc == -1) { csrc = mypcol; };


 mb = descA[MB_];
 nb = descA[NB_];

 lrindx1 = -1; 
 lrindx2 = -1;

 is_found = 0;  
 for(i=ia; i <= ia + m-1; i++) {
   iproc = Cindxg2p( i, mb,myprow,rsrc,nprow);
   is_found = (iproc == myprow);
   if (is_found) {
     lrindx1 = Cnumroc(i,mb,myprow,rsrc,nprow);
     break;
   };
 };

 is_found = 0;
 for(i=(ia+m-1); i >= ia; i--) {
   iproc = Cindxg2p( i, mb,myprow,rsrc,nprow);
   is_found = (iproc == myprow);
   if (is_found) {
     lrindx2 = Cnumroc(i,mb,myprow,rsrc,nprow);
     break;
   };
 };

 *msize = 0;
 *lr1 = lrindx1;
 *lr2 = lrindx2;
 if ((1 <= lrindx1) && (lrindx1 <= lrindx2) ) {
    *msize = lrindx2 - lrindx1 + 1;
    if (!(lrindx2 <= descA[LLD_])) {
      printf("myprow %d mypcol %d lrindx2 %d  descA[LLD_] %d\n",
              myprow, mypcol, lrindx2, descA[LLD_] );
      printf("ia %d m %d mb %d descA[RSRC_] %d lrindx1 %d \n",
              ia,   m,  descA[MB_], descA[RSRC_], lrindx1 );
    };
    assert( lrindx2 <= descA[LLD_] );
 };

 lcindx1 = -1;
 lcindx2 = -1;

 is_found = 0;
 for(j=ja; j <= ja +n-1; j++) {
   iproc = Cindxg2p(j,nb,mypcol,csrc,npcol);
   is_found = (iproc == mypcol);
   if (is_found) {
     lcindx1 = Cnumroc(j,nb,mypcol,csrc,npcol);
     break;
   };
 };

 is_found = 0;
 for(j=ja+n-1; j >= ja; j--) {
   iproc = Cindxg2p(j,nb,mypcol,csrc,npcol);
   is_found = (iproc == mypcol);
   if (is_found) {
     lcindx2 = Cnumroc(j,nb,mypcol,csrc,npcol);
     break;
   };
 };

 *nsize = 0;
 *lc1 = lcindx1;
 *lc2 = lcindx2;
 if ((1 <= lcindx1) && (lcindx1 <= lcindx2)  ) {
   *nsize = (lcindx2 - lcindx1 + 1);
 };

 if (idebug >= 1) {
   printf("local_extent:msize %d  lr1 %d lr2 %d  nsize %d lc1 %d lc2 %d\n",
       *msize, *lr1, *lr2,    *nsize,*lc1,*lc2 );
   printf("lrindx1 %d lrindx2 %d lcindx1 %d lcindx2 %d \n",
           lrindx1,   lrindx2,   lcindx1,   lcindx2 );
 };

 return;
}


#ifdef __cplusplus
extern "C" 
#endif
void local_extent(int m,int n,int ia,int ja,int *descA, 
                  int *msize, int *nsize, 
                  int *lr1, int *lc1, int *lr2, int *lc2 )
{
  const int idebug = 0;
  int isok;

  int lr1_info,lr2_info, lc1_info,lc2_info, msize_info,nsize_info;
  int lr1_loop,lr2_loop, lc1_loop,lc2_loop, msize_loop,nsize_loop;
  int mb,nb,rsrc,csrc;
  int icontxt, nprow,npcol,  myprow,mypcol;

  const int use_local_extent_loop = FALSE;

  if (use_local_extent_loop) {
    local_extent_loop(m,n,ia,ja,descA,  
           msize,nsize,   lr1,lc1,lr2,lc2);
    return;
  };



  icontxt = descA[CTXT_];
  Cblacs_gridinfo( icontxt, &nprow,&npcol, &myprow,&mypcol);

  mb = descA[MB_];
  nb = descA[NB_];
  rsrc = descA[RSRC_];
  csrc = descA[CSRC_];


  local_extent_info( m,n, ia,ja, descA, 
      &msize_info, &nsize_info,
      &lr1_info, &lc1_info, &lr2_info, &lc2_info );

  if (idebug >= 1) {
    local_extent_loop( m,n,  ia,ja, descA,
          &msize_loop, &nsize_loop,
          &lr1_loop, &lc1_loop,   &lr2_loop, &lc2_loop );

    isok = (msize_loop == msize_info) &&
           (nsize_loop == nsize_info);


    if (isok) {
      if (msize_loop >= 1) {
           isok = isok && 
             (lr1_loop == lr1_info) &&
             (lr2_loop == lr2_info); 
      };
      if (nsize_loop >= 1) {
           isok = isok &&
             (lc1_loop == lc1_info) &&
             (lc2_loop == lc2_info) ;
      };
    };
    if (!isok) {
      printf("local_extend:\n");
      printf("myprow %d mypcol %d nprow %d npcol %d \n",
              myprow,   mypcol,   nprow,   npcol );
      printf("m %d n %d ia %d ja %d \n",
              m,   n,   ia,   ja );
      printf("mb %d nb %d rsrc %d csrc %d \n",
              mb,   nb,   rsrc,   csrc );

      printf("msize_loop %d msize_info %d \n", 
              msize_loop,   msize_info );
      printf("nsize_loop %d nsize_info %d \n", 
              nsize_loop,   nsize_info );
      printf("lr1_loop %d lr1_info %d \n", 
              lr1_loop,   lr1_info );
      printf("lr2_loop %d lr2_info %d \n", 
              lr2_loop,   lr2_info );
      printf("lc1_loop %d lc1_info %d \n", 
              lc1_loop,   lc1_info );
      printf("lc2_loop %d lc2_info %d \n", 
              lc2_loop,   lc2_info );
    };


  };

  *lr1 = lr1_info;
  *lr2 = lr2_info;
  *lc1 = lc1_info;
  *lc2 = lc2_info;

  *msize = msize_info;
  *nsize = nsize_info;


  return;
}



