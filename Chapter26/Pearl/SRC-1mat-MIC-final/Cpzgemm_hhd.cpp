#include "ooclu.h"

#define dC(i,j)  (((cuDoubleComplex *)C) + IDX2F((i),(j),descC[LLD_]))
#define dA(i,j)  (((cuDoubleComplex *)dAtmp) + IDX2F((i),(j),descAtmp[LLD_]))
#define dB(i,j)  (((cuDoubleComplex *)dBtmp) + IDX2F((i),(j),descBtmp[LLD_]))

static int nbytesAtmp = 0;
static double *Atmp = 0;
static cuDoubleComplex *dAtmp = 0;

static int nbytesBtmp = 0;
static double *Btmp = 0;
static cuDoubleComplex *dBtmp = 0;


#ifdef __cplusplus
extern "C"
#endif
void Cpzgemm_hhd(char transA, char transB, int m, int n, int k,
        double *alpha,  double *A, int ia,int ja,int *descA,
	                double *B, int ib,int jb,int *descB,
        double *beta,   cuDoubleComplex *C, int ic,int jc,int *descC )
{

/*
 * simulate PZGEMM but A, B are distributed matrices on host
 * C is distributed memory on GPU
 *
 * C <- beta * C + alpha * op(A)  * op(B)
 */


  const int use_MallocHost = FALSE;
  const int use_broadcast = FALSE;
  const int idebug = 0;

  int ip,isok;

  int Locp, Locq, lld, info;
  int mb,nb, ic_proc, jc_proc;
  int icontxt, nprow,npcol,myprow,mypcol;
  int k1, k2,   mm,nn,kk;
  int notransA, notransB;
  int iia,jja,  iib,jjb, iic,jjc;
  int is_k_small, has_work;
  int k_max, inc1, inc2;

  int elemSize = sizeof( cuDoubleComplex );
  cublasStatus  cu_status;
  size_t nbytesA, nbytesB;

  int lrA1,lcA1, lrA2,lcA2;
  int lrB1,lcB1, lrB2,lcB2;
  int lrC1,lcC1, lrC2,lcC2;
  int LocpA, LocqA,  LocpB, LocqB, LocpC, LocqC;

  cuDoubleComplex zalpha, zbeta;
  int ldA,ldB,ldC;
  int ldAtmp, ldBtmp;
  int isizeAtmp, isizeBtmp;

  int descBtmp[DLEN_];
  int descAtmp[DLEN_];
  int rsrc, csrc;
  int mbnb;

  double one_[REAL_PART+IMAG_PART+1];
  double zero_[REAL_PART+IMAG_PART+1];
  double *one = &(one_[0]);
  double *zero = &(zero_[0]);

  double new_beta_[REAL_PART+IMAG_PART+1];
  double *new_beta = &(new_beta_[0]);
  double new_alpha_[REAL_PART+IMAG_PART+1];
  double *new_alpha = &(new_alpha_[0]);


  
  one[REAL_PART] = 1.0;
  one[IMAG_PART] = 0.0;

  zero[REAL_PART] = 0.0;
  zero[IMAG_PART] = 0.0;


  has_work = (m >= 1) && (n >= 1) && (k >= 1);
  if (!has_work) {
      return;
      };



 PROFSTART("Cpzgemm_hhd");


  icontxt = descC[CTXT_];
  Cblacs_gridinfo( icontxt, &nprow,&npcol,&myprow,&mypcol);
  assert( nprow >= 1);
  assert( npcol >= 1);
  assert( (0 <= myprow) && (myprow < nprow));
  assert( (0 <= mypcol) && (mypcol < npcol));


  notransA = (transA == 'N') || (transA == 'n');
  notransB = (transB == 'N') || (transB == 'n');

  mbnb = MIN( descC[MB_], descC[NB_] );
  k_max = mbnb;
  is_k_small = (k <= k_max );

  if (is_k_small) {


     /* 
      * broadcast copy into buffer on GPU
      * perform GEMM on GPU
      */


    /*
     * setup temp buffers
     * Atmp(m,:) should be aligned to C(ic:(ic+m-1),:)
     * Btmp(:,n) should be aligned to C(:, jc:(jc+n-1))
     */
    mb = descC[MB_];
    nb = descC[NB_];
    ic_proc = Cindxg2p( ic, mb, myprow, descC[RSRC_], nprow);
    jc_proc = Cindxg2p( jc, nb, mypcol, descC[CSRC_], npcol);


    rsrc = ic_proc; 
    csrc = mypcol;

    LocpA = Cnumroc( m, mb, 0,0,nprow );
    LocqA = Cnumroc( k, nb, 0,0,npcol );
    lld = MAX(1,LocpA);
    isizeAtmp = lld *  k;

    info = 0;
    Cdescinit( descAtmp, m,k,  
           mb, k,   rsrc, csrc, icontxt, lld, &info );
    assert( info == 0);


    rsrc = myprow;
    csrc = jc_proc;


    LocpB = Cnumroc( k, mb, 0,0,nprow);
    LocqB = Cnumroc( n, nb, 0,0,npcol);
    lld = MAX(1,k);
    isizeBtmp = lld * MAX(1,LocqB);

    info = 0;
    Cdescinit( descBtmp, k,n,
               k, nb,   rsrc, csrc, icontxt, lld, &info);
    assert( info == 0);


     nbytesA = elemSize;
     nbytesA *= isizeAtmp;

     if (nbytesAtmp <  nbytesA) {

       /*
        * Current buffer not large enough
        * Free and reallocate
        */
      if (Atmp != 0) {
         if (use_MallocHost) {
             FreeHost(Atmp); 
         }
         else {
           free(Atmp);
         };

         Atmp = 0;
         };

       if (dAtmp != 0) {
         CUBLAS_FREE( dAtmp );
         dAtmp = 0;
       };

         

     if (use_MallocHost) {
       Atmp = (double *) MallocHost( nbytesA );
     }
     else {
       Atmp = (double *) malloc( nbytesA );
     };
     assert( Atmp != 0 );

     CUBLAS_MALLOC( dAtmp, isizeAtmp, elemSize );
     assert( dAtmp != 0);

     nbytesAtmp = nbytesA;
     };




     nbytesB = elemSize;
     nbytesB *= isizeBtmp;


     if (nbytesBtmp < nbytesB ) {
       /*
        * free and reallocate larger buffer
        */

      if (Btmp != 0) {
         if (use_MallocHost) {
             FreeHost(Btmp); 
         }
         else {
           free(Btmp);
         };

         Btmp = 0;
         };

       CUBLAS_FREE( dBtmp );
       dBtmp = 0;



     if (use_MallocHost) {
       Btmp = (double *) MallocHost( nbytesB );
     }
     else {
       Btmp = (double *) malloc( nbytesB );
     };
     assert( Btmp != 0 );

     CUBLAS_MALLOC( dBtmp, isizeBtmp, elemSize );
     assert( dBtmp != 0);

     nbytesBtmp = nbytesB;
     };





    local_extent( m,n, ic,jc,descC,   
                   &LocpC, &LocqC, &lrC1,&lcC1, &lrC2,&lcC2 );

     /*
      * perform broadcast copy
      */

     /*
      * set descAtmp[CSRC_] = -1, descBtmp[RSRC_] = -1
      * as replicated data  for broadcast copy by PBLAS
      */
    PROFSTART("gemm:copy");
    if (use_broadcast) {
      char scope = 'R';
      char top = ' ';

      descAtmp[CSRC_] = jc_proc;
      iia = 1; 
      jja = 1;
      new_alpha = one;
      new_beta = zero;
      scalapack_pzgeadd( &transA, &m, &k, 
          new_alpha, A, &ia,&ja,descA,
          new_beta,  Atmp, &iia,&jja,descAtmp );

      Locp = isizeAtmp;
      Locq = 1;
      lld = Locp;
      scope = 'R';
      if (mypcol == jc_proc) {
        scalapack_zgebs2d(&icontxt, &scope, &top, 
            &Locp, &Locq, Atmp, &lld );
        }
      else {
        scalapack_zgebr2d(&icontxt,&scope,&top,
            &Locp,&Locq, Atmp, &lld,
               &myprow, &jc_proc );
        };

      descAtmp[CSRC_] = mypcol;
     /*
      * copy buffer to GPU
      */

     inc1 = 1;
     inc2 = 1;
     cu_status = cublasSetVector(isizeAtmp, elemSize,
                     Atmp, inc1, dAtmp, inc2 );
     CHKERR( cu_status );


      iib = 1;
      jjb = 1;
      descBtmp[RSRC_] = ic_proc;

      scalapack_pzgeadd( &transB, &k, &n,
            new_alpha, B, &ib, &jb, descB,
            new_beta,  Btmp, &iib, &jjb, descBtmp );

      Locp = isizeBtmp;
      Locq = 1;
      lld = Locp;
      scope = 'C';
      if (myprow == ic_proc) {
        scalapack_zgebs2d(&icontxt, &scope, &top,
            &Locp, &Locq, Btmp, &lld );
        }
      else {
        scalapack_zgebr2d(&icontxt, &scope, &top,
            &Locp, &Locq, Btmp, &lld,
                &ic_proc, &mypcol );
      };

      descBtmp[RSRC_] = myprow;
     /*
      * copy buffer to GPU
      */

     inc1 = 1;
     inc2 = 1;
     cu_status = cublasSetVector(isizeBtmp, elemSize,
                     Btmp, inc1, dBtmp, inc2 );
     CHKERR( cu_status );
    }
    else {

     descAtmp[CSRC_] = -1;
     descBtmp[RSRC_] = -1;

     iia = 1; 
     jja = 1;
     new_alpha = one;
     new_beta = zero;

     scalapack_pzgeadd( &transA, &m, &k, 
             new_alpha, A,&ia,&ja,descA,
             new_beta, Atmp,&iia,&jja,descAtmp );
     /*
      * copy buffer to GPU
      */

     inc1 = 1;
     inc2 = 1;
#ifdef USE_CUBLASV2
     {
       cublasStatus_t ierr;
       ierr = cublasSetVectorAsync(isizeAtmp, elemSize,
                     Atmp, inc1, dAtmp, inc2, cublas_get_stream() );
       assert( ierr == CUBLAS_STATUS_SUCCESS );
     }
#else
     cu_status = cublasSetVector(isizeAtmp, elemSize,
                     Atmp, inc1, dAtmp, inc2 );
     CHKERR( cu_status );
#endif

     iib = 1; jjb = 1;
     new_alpha = one;
     new_beta = zero;
     scalapack_pzgeadd( &transB, &k, &n,
            new_alpha,  B, &ib, &jb, descB,
            new_beta, Btmp, &iib, &jjb, descBtmp );
     /*
      * copy buffer to GPU
      */

     inc1 = 1;
     inc2 = 1;
#ifdef USE_CUBLASV2
     {
       cublasStatus_t ierr;
       ierr = cublasSetVectorAsync( isizeBtmp, elemSize,
               Btmp, inc1, dBtmp, inc2, cublas_get_stream() );
       assert( ierr == CUBLAS_STATUS_SUCCESS );
     }
#else
     cu_status = cublasSetVector(isizeBtmp, elemSize,
                     Btmp, inc1, dBtmp, inc2 );
     CHKERR( cu_status );
#endif
     
    };

    PROFEND("gemm:copy");
     /*
      * Treat as replicated storage
      */
     descAtmp[CSRC_] = mypcol;
     descBtmp[RSRC_] = myprow;



     /*
      * check extend of Atmp
      */

     iia = 1; 
     jja = 1;
     local_extent( m,k, iia,jja,descAtmp, 
                  &LocpA, &LocqA, 
                  &lrA1,&lcA1, &lrA2,&lcA2 );




     /*
      * check extend of Btmp
      */

     iib = 1;
     jjb = 1;
     local_extent( k,n, iib,jjb,descBtmp, 
                   &LocpB, &LocqB, 
                   &lrB1,&lcB1, &lrB2,&lcB2 );
     /*
      * Perform local GEMM operation on GPU
      */



   

     mm = LocpC;
     nn = LocqC;;
     kk = k;

     has_work = (mm >= 1) && (nn >= 1) && (kk >= 1);
     if (has_work) {
        assert( mm == LocpA );
        assert( kk == LocqA );

        assert( kk == LocpB );
        assert( nn == LocqB );



     ldAtmp = descAtmp[LLD_];
     ldBtmp = descBtmp[LLD_];
     ldA = ldAtmp;
     ldB = ldBtmp;
     ldC = descC[LLD_];


     zalpha = make_cuDoubleComplex( alpha[REAL_PART], alpha[IMAG_PART]);
     zbeta  = make_cuDoubleComplex( beta[REAL_PART], beta[IMAG_PART]);

     if (idebug >= 2) {
       printf("Cpzgemm_hhd: mm %d nn %d kk %d \n", mm,nn,kk);
       printf("ic %d jc %d lrC1 %d lcC1 %d lrC2 %d lcC2 %d\n",
               ic,   jc,   lrC1,   lcC1,   lrC2,   lcC2 );
       printf("descC: M %d N %d MB %d NB %d LLD %d \n",
               descC[M_], descC[N_], descC[MB_], descC[NB_],
               descC[LLD_] );
       printf("lrA1 %d lcA1 %d lrA2 %d lcA2 %d descAtmp[LLD_] %d \n",
               lrA1,   lcA1,   lrA2,   lcA2,   descAtmp[LLD_] );
       printf("lrB1 %d lcB1 %d lrB2 %d lcB2 %d descBtmp[LLD_] %d \n",
               lrB1,   lcB1,   lrB2,   lcB2,   descBtmp[LLD_] );
       printf("isizeAtmp %d isizeBtmp %d \n",
               isizeAtmp,   isizeBtmp );
     };




#ifdef USE_FAKE_CUBLAS
/*
 * debug
 */
     if (idebug >= 1) {
          cuDoubleComplex *pAtmp = (cuDoubleComplex *) Atmp;
          cuDoubleComplex *pBtmp = (cuDoubleComplex *) Btmp;
          cuDoubleComplex aij,bij;
          int lld;
          int lrindx,lcindx,rsrc,csrc;

          descAtmp[CSRC_] = (npcol-1);
          lld = descAtmp[LLD_];
          for(int i=1; i <= descAtmp[M_]; i++) {
            for(int j=1; j <= descAtmp[N_];j++) {
                Cinfog2l(i,j,descAtmp,nprow,npcol,myprow,mypcol,
                     &lrindx,&lcindx, &rsrc,&csrc );
                if ((rsrc == myprow) && (csrc == mypcol)) {
                  aij = pAtmp[ IDX2F(lrindx,lcindx,lld) ];
                  printf("Atmp(%d,%d) = (%lf,%lf)\n",
                      i,j, creal(aij), cimag(aij) );
                };
            };
          };



          descBtmp[RSRC_] = (nprow-1);
          lld = descBtmp[LLD_];
          for(int i=1; i <= descBtmp[M_]; i++) {
            for(int j=1; j <= descBtmp[N_];j++) {
                Cinfog2l(i,j,descBtmp,nprow,npcol,myprow,mypcol,
                     &lrindx,&lcindx, &rsrc,&csrc );
                if ((rsrc == myprow) && (csrc == mypcol)) {
                  bij = pBtmp[ IDX2F(lrindx,lcindx,lld) ];
                  printf("Btmp(%d,%d) = (%lf,%lf)\n",
                      i,j, creal(bij), cimag(bij) );
                };
            };
          };


      };
#endif
       CUBLAS_ZGEMM(
           CUBLAS_OP_N,CUBLAS_OP_N, mm,nn,kk,
          zalpha,  (cuDoubleComplex *) dA(lrA1,lcA1), ldAtmp,
                   (cuDoubleComplex *) dB(lrB1,lcB1), ldBtmp,
          zbeta,   (cuDoubleComplex *) dC(lrC1,lcC1), ldC );


       };



     }
  else {
     /*
      *  split operation and use recursion
      */

    k1 = k_max;
    k2 = k - k1;

    Cpzgemm_hhd( transA, transB, m,n,k1,
        alpha,   A, ia,ja,descA,
                 B, ib,jb,descB,
        beta,    C, ic,jc,descC );




    if (notransA) {
      iia = ia;  jja = ja + k1;
      }
    else {
      jja = ja; iia = ia + k1;
      };

    if (notransB) {
      iib = ib + k1; jjb = jb;
      }
    else {
      jjb = jb + k1; iib = ib;
      };

    new_beta = one;
    Cpzgemm_hhd( transA, transB, m,n,k2,
           alpha, A, iia,jja,descA,
                  B, iib,jjb,descB,
           new_beta,  C, ic,jc,  descC );

    };


 PROFEND("Cpzgemm_hhd");

}




#ifdef __cplusplus
extern "C"
#endif
void pzgemm_hhd( char *transA, char *transB,
          int *m, int *n, int *k,
          double *alpha, double *A, int *ia, int *ja, int *descA,
                         double *B, int *ib, int *jb, int *descB,
          double *beta,  cuDoubleComplex *C, int *ic, int *jc, int *descC )
{

Cpzgemm_hhd( *transA, *transB,
       *m, *n, *k,
       alpha, A, *ia, *ja, descA,
              B, *ib, *jb, descB,
       beta,  C, *ic, *jc, descC );
}




