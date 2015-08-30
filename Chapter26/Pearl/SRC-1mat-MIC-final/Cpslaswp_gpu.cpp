#include "ooclu.h"

#define ipiv(i)  ipiv_[IDX1F(i)]

#ifdef __cplusplus
extern "C"
#endif
void
Cpslaswp_gpu( char direc, char rowcol, int n, 
              float *A, int ia, int ja, int *descA, 
              int k1, int k2, int *ipiv_ )
{
  /*
   * Note vector ipiv(:) is tied to the distribution of A
   * dimension ipiv is   Locr(M_A) + MB_A
   */

  int nprow,npcol,myprow,mypcol;
  int ip,jp;
  int i,j, iia,jja, icurrow,icurcol;

  int is_forward;
  int is_row;

  if (n <= 0) {
    return;
  };

  is_row = (rowcol == 'R') || (rowcol == 'r');
  is_forward =  (direc == 'F') || (direc == 'f');

  Cblacs_gridinfo( descA[CTXT_], &nprow,&npcol, &myprow,&mypcol );
  assert( nprow >= 1);
  assert( npcol >= 1);
  assert( (0 <= myprow) && (myprow < nprow));
  assert( (0 <= mypcol) && (mypcol < npcol));
/*
!      IF( LSAME( ROWCOL, 'R' ) ) THEN
!         IF( LSAME( DIREC, 'F' ) ) THEN
!            CALL INFOG2L( K1, JA, DESCA, NPROW, NPCOL, myprow, mypcol,
!     $                    IIA, JJA, ICURROW, ICURCOL )
!            DO 10 I = K1, K2
!               IP = IPIV( IIA+I-K1 )
!               IF( IP.NE.I )
!     $            CALL PZSWAP( N, A, I, JA, DESCA, DESCA( M_ ), A, IP,
!     $                         JA, DESCA, DESCA( M_ ) )
!   10       CONTINUE
*/
      if (is_row && is_forward ) {
             Cinfog2l(k1,ja,descA,nprow,npcol,myprow,mypcol,
                 &iia,&jja,   &icurrow,&icurcol );
             for(i=k1; i <= k2; i++) {
               ip = ipiv( iia + i - k1 );
               if (ip != i) {
                  Cpsswap_gpu(n, A, ia,ja,descA, descA[M_], A, ip,
                               ja, descA, descA[M_] );
               };
             }; /* end for i */
         };
/*
!         ELSE
!            CALL INFOG2L( K2, JA, DESCA, NPROW, NPCOL, myprow, mypcol,
!     $                    IIA, JJA, ICURROW, ICURCOL )
!            DO 20 I = K2, K1, -1
!               IP = IPIV( IIA+I-K1 )
!               IF( IP.NE.I )
!     $            CALL PZSWAP( N, A, I, JA, DESCA, DESCA( M_ ), A, IP,
!     $                         JA, DESCA, DESCA( M_ ) )
!   20       CONTINUE
!         END IF
*/
         if (is_row && (!is_forward)) {
           Cinfog2l(k2,ja,descA,nprow,npcol,myprow,mypcol,
                       &iia,&jja,  &icurrow,&icurcol );
           for(i=k2; i >= k1; i-- ) {
             ip = ipiv( iia + i - k1 );
             if (ip != i) {
               Cpsswap_gpu( n,A,i,ja,descA, descA[M_],  A, ip,
                            ja, descA, descA[M_] );
             };
           }; /* end for i */
         };
/*
!      ELSE
!         IF( LSAME( DIREC, 'F' ) ) THEN
!            CALL INFOG2L( IA, K1, DESCA, NPROW, NPCOL, myprow, mypcol,
!     $                    IIA, JJA, ICURROW, ICURCOL )
!            DO 30 J = K1, K2
!               JP = IPIV( JJA+J-K1 )
!               IF( JP.NE.J )
!     $            CALL PZSWAP( N, A, IA, J, DESCA, 1, A, IA, JP,
!     $                         DESCA, 1 )
!   30       CONTINUE
*/
         if ( (!is_row) && is_forward) {
               Cinfog2l(ia,k1,descA,  nprow,npcol,myprow,mypcol,
                           &iia, &jja,  &icurrow, &icurcol );
               for(j=k1; j <= k2; j++ ) {
                 jp = ipiv( jja+j-k1 );
                 if (jp != j) {
                   Cpsswap_gpu( n,A,ia,j,descA,1,  A,ia,jp,
                                descA, 1 );
                 };
               }; /* end for j */
           };
/*
!         ELSE
!            CALL INFOG2L( IA, K2, DESCA, NPROW, NPCOL, myprow, mypcol,
!     $                    IIA, JJA, ICURROW, ICURCOL )
!            DO 40 J = K2, K1, -1
!               JP = IPIV( JJA+J-K1 )
!               IF( JP.NE.J )
!     $            CALL PZSWAP( N, A, IA, J, DESCA, 1, A, IA, JP,
!     $                         DESCA, 1 )
!   40       CONTINUE
!         END IF
*/
           if ( (!is_row) && (!is_forward))  {
             Cinfog2l( ia,k2,descA, nprow,npcol,myprow,mypcol,
                          &iia, &jja, &icurrow, &icurcol );
             for( j=k2; j >= k1; j--) {
               jp = ipiv( jja + j-k1);
               if (jp != j) {
                 Cpsswap_gpu( n,A,ia,j,descA, 1, A, ia,jp, 
                              descA, 1 );
               };
             }; /* end for j */
           };
/*
!      END IF
*/
         return;
}
