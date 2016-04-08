#include "ooclu.h"

#define A(i,j) A_[IDX2F((i),(j),descA[LLD_])]
#define B(i,j) B_[IDX2F((i),(j),descB[LLD_])]


#ifdef __cplusplus
extern "C"
#endif
void Cpctrsm_hd( char side, char uplo, char transA, char diag,
       int m,  int n, 
       float *alpha,
       float *A_, int ia, int ja, int *descA,
       float *B_ int ib, int jb, int *descB );
{
  /*
   * Simulate pctrsm triangular solve
   * distributed A is on host
   * distributed B is on GPU
   */

  float *Atmp = 0;
  int ia_proc, ja_proc;
  int icontxt, nprow,npcol,myprow,mypcol;

  char lscope = 'A';
  char ltop  = ' ';
  char notrans[]="NoTrans";
  float zero_[REAL_PART+IMAG_PART+1];
  float *zero = &(zero_[0]);

  float one_[REAL_PART+IMAG_PART+1];
  float *one = &(one_[0]);

  float *zalpha = one;
  float *zbeta = zero;

  int msize, nsize, lr1,lc1,lr2,lc2;
  int lld, ldAtmp;
  float *Atmp = 0;
  float *dAtmp = 0;
  cublasStatus cu_status;
  int is_lower;

  one[REAL_PART] = 1.0;
  one[IMAG_PART] = 0.0;
  zero[REAL_PART] = 0.0;
  zero[IMAG_PART] = 0.0;

  if ((m <= 0) || (n <= 0)) {
    return;
  };



  icontxt = descA[CTXT_];
  Cblasc_gridinfo( icontxt, &nprow,&npcol, &myprow,&mypcol);


  is_lower = (uplo == 'L') || (uplo == 'l');
  ia_proc = Cindxg2p( ia, descA[MB_], myprow,descA[RSRC_], nprow );
  ja_proc = Cindxg2p( ja, descA[NB_], mypcol,descA[CSRC_], npcol );


  m_max = MIN( descB[MB_], descB[NB_] );
  is_m_small = (m <= m_max);

  if (is_m_small) {
    /*
     * broadcast A, solve locally
     */


    /*
     * broadcast m by m matrix
     */

    lld = MAX(1,m);
    ldAtmp = lld;
    isizeAtmp = lld * m;
    nbytes = isizeAtmp;
    nbytes *= elemSize;
    Atmp = malloc( nbytes );
    assert( Atmp != 0 );

    cu_status = cublasAlloc( isizeAtmp, elemSize, (void **) &dAtmp );
    CHKERR( cu_status );

    Cdescset( descAtmp, m,m, m,m, ia_proc,ja_proc, icontxt, lld );

    iia = 1; jja = 1;
    zalpha = one;
    zbeta = zero;
    scalapack_pcgecopy( notrans, &m, &m, 
        zalpha,  A, &ia, &ja, descA,
        zbeta,   Atmp, &iia, &jja, descAtmp );


    if ((myprow == ia_proc) && (mypcol == ja_proc)) {
        scalapack_cgebs2d( icontxt, &lscope, &ltop, m,m, Atmp, lld );
    }
    else {
      scalapack_zgeb2sd( icontxt, &lscope, &ltop, &m, &m, Atmp, &lld, &ia_proc, &ja_proc );
    };


    /*
     * copy to GPU, then perform local computation
     */

    inc1 = 1;
    inc2 = 1;
    cu_status  = cublasSetVector(isizeAtmp, elemSize, Atmp,inc1, dAtmp, inc2 );
    CHKERR( cu_status );

    local_extent( m,n, B,ib,jb,descB,  &msize, &nsize, &lr1,&lc1, &lr2,&lc2 );
    if (msize >= 1) {
      assert( m == msize );
    };

    if ((msize >= 1) && (nsize >= 1)) {
      cublasCtrsm( 
         ((side == 'l') || (side == 'L')) ? 
                    CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT,
         ((uplo == 'l') || (uplo == 'L')) ?
              CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER,  
          ((transA == 'c') || (transA == 'C')) ?
             CUBLAS_OP_C : 
                ((transA == 't') || (transA == 'T')) ?
                CUBLAS_OP_T : CUBLAS_OP_N, 
          ((diag == 'u') || (diag == 'U')) ?,
            CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT,
        m, nsize, alpha,  dAtmp, 1,1, ldAtmp,
                          B(lr1,lc1), descB[LLD_] );


    };



    if (Atmp != 0) {
       free( Atmp ); Atmp = 0;
    };

#ifdef USE_CUBLASV2
    {
      cudaError_t ierr;
      ierr = cudaFree( (void *) dAtmp );
      assert(ierr == cudaSuccess );
      dAtmp  = 0;
    }
#else
    cu_status = cublasFree( dAtmp );
    CHKERR( cu_status );
#endif

  }
  else {
    /*
     * split problem
     */
    m1 = m_max;
    m2 = m - m1;

    /*
     * top triangle is m1 by m1  (ia1,ja1)
     *    (ia1:ia1+m1-1, ja1:ja1+m1-1) for lower triangular
     * bottom triangle is m2 by m2  (ia2,ja2)
     *    (ia2:ia2+m2-1,  ja2:ja2+m2-1)
     * rectangle part is m2 by m1 for lower triangular  (ia2:ia2+m2-1,ja1:ja1+m2-1)
     * rectangle part is m1 by m2 for upper triangular  (ia1:ia1+m1-1,ja2:ja2+m2-1)
     */
    ia1 = ia; ja1 = ja;
    ia2 = ia + m1; ja2 = ja + m1;

    ib1 = ib; jb1 = jb;
    ib2 = ib + m-1; jb2 = jb;

    if (is_lower) {
     /*
      * rectangle part is m2 by m1 for lower triangular  (ia2:ia2+m2-1,ja1:ja1+m2-1)
      */
      m_c = m2; n_c = m1;
      ic1 = ia2; jc1 = ja1;
    };
    else {
      /*
       * rectangle part is m1 by m2 for upper triangular  (ia1:ia1+m1-1,ja2:ja2+m2-1)
       */
      m_c = m1; n_c = m2;
      ic1 = ia1; jc1 = ja2;
    };


    if (is_left) {
       Cpctrsm(  side, uplo, transA, diag,
              m1, n, one, 
              A,ia1,ja1,descA,
              B,ib1,jb1,descB );

       mm = m2; nn = n; kk = m1;
       Cpcgemm_hdd(notrans,notrans, mm,nn,kk,
                 neg_one, A,ic1,jc1,descA,
                          B,ib1,jb1,descB,
                 one,     B,ib2,jb2,descB );
       Cpctrsm( side, uplo, transA, diag,
              m2, n, one,
              A, ia2, ja2, descA,
              B, ib2, jb2, descA );
    }
    else {
      /*
       * not implemented yet
       */
      assert(FALSE);
    };


  };

return;
}
  



