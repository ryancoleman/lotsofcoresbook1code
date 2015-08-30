#include "ooclu.h"

#define dC(i,j)  (((double *)C) + IDX2F((i),(j),descC[LLD_]))
#define dA(i,j)  (((double *)dAtmp) + IDX2F((i),(j),descAtmp[LLD_]))
#define dB(i,j)  (((double *)dBtmp) + IDX2F((i),(j),descBtmp[LLD_]))

static int nbytesAtmp = 0;
static double *Atmp = 0;
static double *dAtmp = 0;

static int nbytesBtmp = 0;
static double *Btmp = 0;
static double *dBtmp = 0;



/*
 * Recursive helper routine
 * arrays dAtmp, dBtmp, C are all arrays already in GPU device memory
 */
#ifdef __cplusplus
extern "C"
#endif
void Cpdsyrk_helper( char uplo, int n, int k,
                     double *alpha,  double *dAtmp,  int ia, int ja, int *descAtmp,
                                     double *dBtmp,  int ib, int jb, int *descBtmp,
                     double *beta,   double *C,      int ic, int jc, int *descC,
                     int myprow, int mypcol, int nprow, int npcol )
{

double zalpha = (double) alpha[REAL_PART];
double zbeta = (double) beta[REAL_PART];

int lrC1, lcC1, lrC2, lcC2;
int lrA1, lcA1, lrA2, lcA2;
int lrB1, lcB1, lrB2, lcB2;
int LocpA, LocqA;
int LocpB, LocqB;
int LocpC, LocqC;
int n_small;

int ic11,jc11,  ic12,jc12, ic21,jc21, ic22,jc22;
int ia11,ja11,  ia21,ja21;
int ib11,jb11,  ib12,jb12;


int is_lower = (uplo == 'L') || (uplo == 'l');

int ldAtmp = descAtmp[LLD_];
int ldBtmp = descBtmp[LLD_];
int ldC    = descC[LLD_];

int mb = descC[MB_];
int nb = descC[NB_];
int mm, nn, kk;

int n1, n2, jc_end;

int iic, jjc, rsrc1, csrc1, rsrc2, csrc2;
int is_same_block;
int is_ok, is_mine;
int jstart, jend, jsize, jcfirst, jfinal;
int istart, iend, isize, icfirst, ifinal;
int iia,jja, iib,jjb;

int is_diagonal;
int has_work;


      /*
       * check whether the problem is so small that it is within
       * a single block
       */

     ifinal = ic + n - 1;
     jfinal = jc + n - 1;

     n_small = MAX( nprow*mb, npcol*nb);

     if (n <= n_small) {

     if (is_lower) {
      jcfirst = Ciafirst(jc, nb, mypcol, descC[CSRC_], npcol);
      for(jstart=jcfirst; jstart <= jfinal; jstart = jend + (npcol-1)*nb + 1) {
           /* 
            * set jend on end of block boundary
            */
           assert( mypcol == 
              Cindxg2p(jstart, nb, mypcol, descC[CSRC_],npcol)  );

           jend = jstart + (nb - MOD(jstart-1,nb) - 1);
           jend = MIN(jend, jfinal);
           jsize = jend - jstart + 1;

           assert( mypcol == 
              Cindxg2p(jend, nb, mypcol, descC[CSRC_],npcol)  );

           /*
            * find 1st global row index in lower triangular part
            * that is on this local processor
            */
         
           iic = ic + (jstart-jc);
           icfirst = Ciafirst(iic, mb, myprow, descC[RSRC_],nprow);
           istart = icfirst;


           assert( myprow ==
              Cindxg2p(istart,mb,myprow,descC[RSRC_],nprow) );

           is_diagonal = (istart - ic) == (jstart - jc);
           if (is_diagonal) {

              mm = jsize;
              nn = jsize;
              kk = k;

              iia = ia + (istart-ic);
              jja = ja;
              local_extent( mm, kk,  iia, jja, descAtmp,
                  &LocpA, &LocqA, &lrA1, &lcA1,  &lrA2, &lcA2 );

              
              iic = istart;
              jjc = jstart;
              local_extent( mm, nn, iic, jjc, descC,
                  &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );



#ifdef USE_CUBLASV2
              CUBLAS_DSYRK( 
                    (uplo == 'L') || (uplo == 'l') ? 
                               CUBLAS_FILL_MODE_LOWER :
                               CUBLAS_FILL_MODE_UPPER,
                    CUBLAS_OP_N,
                    nn, kk,
                    zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                    zbeta,   (double *) dC(lrC1,lcC1), ldC );
#else
              CUBLAS_DSYRK( uplo, 'N', nn, kk,
                    zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                    zbeta,   (double *) dC(lrC1,lcC1), ldC );
#endif

              /*
               * diagonal block is done,
               * shift to next block
               */
              istart = istart + jsize;
              };

            /*
             * handle remaining off-diagonal part
             */

              mm = ifinal - istart + 1;
              nn = jsize;
              kk = k;

              iic = istart;
              jjc = jstart;
              local_extent( mm, nn, iic, jjc, descC,
                  &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );

              has_work = (LocpC >= 1) && (LocqC >= 1);
              if (has_work) {
                iia = ia + (istart-ic);
                jja = ja;
                local_extent( mm, kk,  iia, jja, descAtmp,
                  &LocpA, &LocqA, &lrA1, &lcA1,  &lrA2, &lcA2 );


                iib = ib;
                jjb = jb + (jstart-jc);
                local_extent(kk, nn, iib,jjb,descBtmp,
                  &LocpB, &LocqB, &lrB1, &lcB1,  &lrB2, &lcB2 );
                
                /*
                 * local sizes
                 */

                mm = LocpC;
                nn = LocqC;
                kk = k;

                assert( LocpC == LocpA );
                assert( LocqC == LocqB );

                assert( LocqA == k );
                assert( LocpB == k );

#ifdef USE_CUBLASV2
                CUBLAS_DGEMM(
                      CUBLAS_OP_N, CUBLAS_OP_N,
                      mm, nn, kk,
                      zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                 (double *) dB(lrB1,lcB1), ldBtmp,
                      zbeta,     (double *) dC(lrC1,lcC1), ldC );
#else
                CUBLAS_DGEMM(
                      'N',  'N',
                      mm, nn, kk,
                            zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                       (double *) dB(lrB1,lcB1), ldBtmp,
                            zbeta,     (double *) dC(lrC1,lcC1), ldC );
#endif
               };
            }; /* end for jstart */

         }
       else {
         /*
          * is upper triangular part
          */
         icfirst = Ciafirst( ic, mb, myprow, descC[RSRC_], nprow );
         ifinal = ic + n - 1;
         jfinal = jc + n - 1;
         for(istart=icfirst; istart <= ifinal; istart = iend + (nprow-1)*mb+1) {
             iend = istart + (mb - MOD(istart-1,mb) - 1 );
             iend = MIN( ifinal, iend );
             isize = iend - istart + 1;

             jjc = jc + (istart-ic);
             jstart = Ciafirst( jjc, nb, mypcol, descC[CSRC_], npcol);

             is_diagonal = (istart - ic) == (jstart - jc);
             if (is_diagonal) {

               
              mm = isize;
              nn = isize;
              kk = k;

              iia = ia + (istart-ic);
              jja = ja;
              local_extent( mm, kk,  iia, jja, descAtmp,
                  &LocpA, &LocqA, &lrA1, &lcA1,  &lrA2, &lcA2 );

              
              iic = istart;
              jjc = jstart;
              local_extent( mm, nn, iic, jjc, descC,
                  &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );



#ifdef USE_CUBLASV2
              CUBLAS_DSYRK( 
                    (uplo == 'L') || (uplo == 'l') ? 
                               CUBLAS_FILL_MODE_LOWER :
                               CUBLAS_FILL_MODE_UPPER,
                    CUBLAS_OP_N,
                    n, k,
                    zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                    zbeta,   (double *) dC(lrC1,lcC1), ldC );
#else
              CUBLAS_DSYRK( uplo, 'N', n, k,
                    zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                    zbeta,   (double *) dC(lrC1,lcC1), ldC );
#endif

              /*
               * diagonal block is done,
               * shift to next block
               */
              jstart = jstart + isize;
              };

            /*
             * handle remaining off-diagonal part
             */

              mm = isize;
              nn = jfinal - jstart + 1;
              kk = k;

              iic = istart;
              jjc = jstart;
              local_extent( mm, nn, iic, jjc, descC,
                  &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );

              has_work = (LocpC >= 1) && (LocqC >= 1);
              if (has_work) {
                iia = ia + (istart-ic);
                jja = ja;
                local_extent( mm, kk,  iia, jja, descAtmp,
                  &LocpA, &LocqA, &lrA1, &lcA1,  &lrA2, &lcA2 );


                iib = ib;
                jjb = jb + (jstart-jc);
                local_extent(kk, nn, iib,jjb,descBtmp,
                  &LocpB, &LocqB, &lrB1, &lcB1,  &lrB2, &lcB2 );
                
                /*
                 * local sizes
                 */

                mm = LocpC;
                nn = LocqC;
                kk = k;

#ifdef USE_CUBLASV2
                CUBLAS_DGEMM(
                      CUBLAS_OP_N, CUBLAS_OP_N,
                      mm, nn, kk,
                      zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                 (double *) dB(lrB1,lcB1), ldBtmp,
                      zbeta,     (double *) dC(lrC1,lcC1), ldC );
#else
                CUBLAS_DGEMM(
                      'N',  'N',
                      mm, nn, kk,
                            zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                       (double *) dB(lrB1,lcB1), ldBtmp,
                            zbeta,     (double *) dC(lrC1,lcC1), ldC );
#endif
               }; /* end if (has_work) */


          }; /* end for istart */
       }; /* end if (is_lower) */
      }
     else  {

        /* 
         * problem still too big
         *
         * split problem into C11 and C21 but try to be aligned 
         * on block boundary
         */
         n1 = (int) (n/2);
         jc_end = jc + n1 -1;
         if (MOD(jc_end,nb) != 0) {
              jc_end = jc_end + (nb - MOD(jc_end-1,nb)-1);
              };

         n1 = jc_end - jc + 1;
         n2 = n - n1;
         is_ok = (n1 >= 1) && (n2 >= 1);
         if (!is_ok) {
              /*
               * something is not quite right
               * fall back setting to something simple
               */
              jc_end = jc + (nb - MOD(jc-1,nb)-1);
              if (MOD(jc,nb) == 0) { jc_end = jc; };
              
              n1 = jc_end - jc + 1;
              n2 = n - n1;
              };
          assert( n1 >= 1 );
          assert( n2 >= 1 );
          assert( n1 + n2 == n );

          /*
           * partition matrix
           *
           * [C11  C12 ]    [A11 ]     
           * [C21  C22 ]    [A21 ]
           *
           * [B11  B12 ]
           *
           */

          ic11 = ic; 
          jc11 = jc;

          ic12 = ic;
          jc12 = jc + n1;

          ic21 = ic + n1; 
          jc21 = jc;

          ic22 = ic + n1;
          jc22 = jc + n1;

          ia11 = ia;
          ja11 = ja;

          ia21 = ia + n1;
          ja21 = ja;


          ib11 = ib;
          jb11 = jb;

          ib12 = ib;
          jb12 = jb + n1;

    
       /*
        * lower triangular part
        *
        * [ C11      ]   <-  beta * [ C11      ]   + alpha * [A11 ]  [B11  B12]
        * [ C21  C22 ]              [ C21  C22 ]             [A21 ]
        *
        * where B11 is a copy of  transpose(A11), B12 is a copy of  transpose(A21)
        *
        *  C11 <--  beta * C11  + alpha * A11 * B11   computed using recursion
        *  C22 <--  beta * C11  + alpha * A21 * B12   computed using recursion
        *
        *  C21 <-- beta * C21 + alpha * A21 * B11  computed using GEMM
        *
        * 
        *
        * upper triangular part
        *
        * [C11  C12 ] <--- beta * [C11  C12]  + alpha * [A11 ]  * [B11  B12]
        * [     C22 ]             [     C22]            [A21 ]
        *
        * where B11 is a copy of transpose(A11), B12 is a copy of transpose(A21)
        *
        *  C11 <---  beta * C11 + alpha * A11 * B11
        *  C22 <---  beta * C22 + alpha * A21 * B12
        *
        *  C12 <---  beta * C12 + alpha + A11 * B12
        *
        */

            zalpha = (double) alpha[REAL_PART];
            zbeta = (double) beta[REAL_PART];





       /*
        * Compute off-diagonal part by GEMM
        */
       if (is_lower) {
 
          /*
           * C21 <-- beta * C21 + alpha * A21 * B11
           *
           * C21 is n2 by n1
           */


          mm = n2;
          nn = n1;
          kk = k;

          local_extent( mm, nn, ic21, jc21, descC,
               &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );

          local_extent( mm, kk, ia21, ja21, descAtmp,
               &LocpA, &LocqA, &lrA1, &lcA1, &lrA2, &lcA2 );

          local_extent( kk, nn, ib11, jb11, descBtmp,
               &LocpB, &LocqB, &lrB1, &lcB1, &lrB2, &lcB2 );

          }
         else {

          /*
           *  C12 <--  beta * C12 + alpha + A11 * B12
           *
           *  C12 is n1 by n2
           */
          mm = n1;
          nn = n2;
          kk = k;

          local_extent( mm, nn, ic12, jc12, descC,
               &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );

          local_extent( mm, kk, ia11, ja11, descAtmp,
               &LocpA, &LocqA, &lrA1, &lcA1, &lrA2, &lcA2 );

          local_extent( kk, nn, ib12, jb12, descBtmp,
               &LocpB, &LocqB, &lrB1, &lcB1, &lrB2, &lcB2 );
          
         };





          /*
           * local piece of C21 for lower part (or C12 for upper part) matrix 
           * on this processor is size LocpC by LocqC
           */
          mm = LocpC;
          nn = LocqC;
          kk = k;

          if ((mm >= 1) && (nn >= 1) && (kk >= 1)) {
                assert( LocpA == LocpC );
                assert( LocqB == LocqC );

#ifdef USE_CUBLASV2
                CUBLAS_DGEMM(
                      CUBLAS_OP_N, CUBLAS_OP_N,
                      mm, nn, kk,
                      zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                 (double *) dB(lrB1,lcB1), ldBtmp,
                      zbeta,     (double *) dC(lrC1,lcC1), ldC );
#else
                CUBLAS_DGEMM(
                      'N',  'N',
                      mm, nn, kk,
                            zalpha,    (double *) dA(lrA1,lcA1), ldAtmp,
                                       (double *) dB(lrB1,lcB1), ldBtmp,
                            zbeta,     (double *) dC(lrC1,lcC1), ldC );
#endif

              };


        /*
         * use recursion to evaluate  C11 and C22
         *  C11 <--  beta * C11  + alpha * A11 * B11   computed using recursion
         *  C22 <--  beta * C11  + alpha * A21 * B12   computed using recursion
         */

             Cpdsyrk_helper( uplo, n1, k,
                       alpha,    dAtmp, ia11, ja11, descAtmp,
                                 dBtmp, ib11, jb11, descBtmp,
                       beta,     C,     ic11, jc11, descC,
                       myprow,mypcol,   nprow,npcol );


             Cpdsyrk_helper( uplo, n2, k,
                       alpha,   dAtmp, ia21, ja21, descAtmp,
                                dBtmp, ib12, jb12, descBtmp,
                       beta,    C,     ic22, jc22, descC,
                       myprow,mypcol,   nprow,npcol );
    };

}


#ifdef __cplusplus
extern "C"
#endif
void Cpdsyrk_hhd(char  	uplo, char  	transA, int m,   int   	n,	int   	k,
                 double *  	alpha, double *  	A, int   	ia, int   	ja, int *  	descA,
                 double *  	beta, double *  	C, int   	ic, int   	jc, int *  	descC )
{
    /*
     * uplo : 'U' or 'L' , standing for upper/lower triangle of C are to be referenced
     *
     * simulate PDSYRK but A is distributed matrices on host
     On entry,  TRANS  specifies the  operation to be performed as
     *          follows:
     *
     *             TRANS = 'N' or 'n'
     *                  sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
     *             TRANS = 'T' or 't'
     *                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ).
     *             TRANS = 'C' or 'c'
     *                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ).
     *
     *
     *     C( ic:(ic+m-1), jc:(jc+n-1)) = beta * C( ic:(ic+m-1), jc:(jc+n-1) ) +   
     *                                    alpha * A( ia:(ia+m-1), ja:(ja+k-1) ) * trans(A( ia:(ia+n-1), ja:(ja+k-1))
     *
     *     C is in GPU device memory
     *     A is distributed in CPU host memory
     */
    char transB='T';

    
    int is_lower = 0;
    int is_upper = 0;


    const int use_MallocHost = FALSE;
    const int use_broadcast = FALSE;
    const int use_copy_Atmp = TRUE;

    const int use_simple = FALSE;
    const int use_recursion = TRUE;

    const int idebug = 1;
    
    int ip,isok;
    
    int Locp, Locq, lld, info;
    int mb,nb, ic_proc, jc_proc;
    int icontxt, nprow,npcol,myprow,mypcol;
    int k1, k2,   mm,nn,kk;
    int notrans;
    int iia,jja,  iib,jjb, iic,jjc;
    int is_k_small, has_work;
    int k_max, inc1, inc2;
    
    int elemSize = sizeof( double );

    
    
    size_t nbytesA, nbytesB;
    
    int ic_diag, jc_diag, ia_diag,ja_diag;

    int lrA1,lcA1, lrA2,lcA2;
    int lrB1,lcB1, lrB2,lcB2;
    int lrC1,lcC1, lrC2,lcC2;
    int LocpA, LocqA,  LocpB, LocqB, LocpC, LocqC;
    
    double zalpha = (double) alpha[REAL_PART]; 
    double zbeta = (double) beta[REAL_PART];
    int ldA,ldB,ldC;
    int ldAtmp, ldBtmp;
    int isizeAtmp, isizeBtmp;
    
    int descBtmp[DLEN_];
    int descAtmp[DLEN_];
    int rsrc, csrc;
    int mbnb;
    int has_diagonal;
    int jstart,jend, jsize, jfinal, jcfirst;
    int nb_A, mb_B;
    
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
    
    
    has_work = (n >= 1)&& (k >= 1);
    if (!has_work) {
        return;
    };


    
    PROFSTART("Cpdsyrk_hhd");
    
    
    icontxt = descC[CTXT_];
    Cblacs_gridinfo( icontxt, &nprow,&npcol,&myprow,&mypcol);
    assert( nprow >= 1);
    assert( npcol >= 1);
    assert( (0 <= myprow) && (myprow < nprow));
    assert( (0 <= mypcol) && (mypcol < npcol));


    is_lower = ( uplo == 'L') || (uplo == 'l');
    is_upper = ( uplo == 'U') || (uplo == 'u');
    assert( is_lower ||  is_upper );
    
    
    notrans = (transA == 'N') || (transA == 'n');
    
    mbnb = MIN( descC[MB_], descC[NB_] );
    k_max = mbnb;
    is_k_small = (k <= k_max );
    
    if (is_k_small) {
        
        /*
         * broadcast copy into buffer on GPU
         * perform SYRK and GEMM on GPU
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
        
        nb_A = k;
        LocpA = Cnumroc( m, mb, 0,0,nprow );
        LocqA = Cnumroc( k, nb_A, 0,0,npcol );
        lld = MAX(1,LocpA);
        isizeAtmp = MAX(1,lld *  k);
        
        info = 0;
        Cdescinit( descAtmp, m,k,
                  mb, nb_A,   rsrc, csrc, icontxt, lld, &info );
        assert( info == 0);
        
        
        rsrc = myprow;
        csrc = jc_proc;
        
        
        mb_B = k;
        LocpB = Cnumroc( k, mb_B, 0,0,nprow);
        LocqB = Cnumroc( n, nb, 0,0,npcol);
        lld = MAX(1,k);
        isizeBtmp = MAX(1,lld * LocqB);
        
        info = 0;
        Cdescinit( descBtmp, k,n,
                  mb_B, nb,   rsrc, csrc, icontxt, lld, &info);
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
            
            if (dBtmp != 0) {
                 CUBLAS_FREE( dBtmp );
                 };
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
        
        
        /*
         * get Atmp and Btmp, broadcast,at most 1 transpose will be done
         */
        
        local_extent( m,n, ic,jc,descC,
                     &LocpC, &LocqC, &lrC1,&lcC1, &lrC2,&lcC2 );

       

#ifndef USE_GPU_DIRECT 
        
        /*
         * set descAtmp[CSRC_] = -1, descBtmp[RSRC_] = -1
         * as replicated data  for broadcast copy by PBLAS
         */
        PROFSTART("syrk:copy");
        if (use_broadcast) {
            char scope = 'R';
            char top = ' ';
            
            descAtmp[CSRC_] = jc_proc;
            iia = 1;
            jja = 1;
            new_alpha = one;
            new_beta = zero;
            scalapack_pdgeadd( &transA, &m, &k,
                              new_alpha, A, &ia,&ja,descA,
                              new_beta,  Atmp, &iia,&jja,descAtmp );
            
            Locp = isizeAtmp;
            Locq = 1;
            lld = Locp;
            scope = 'R';
            if (mypcol == jc_proc) {
                scalapack_dgebs2d(&icontxt, &scope, &top,
                                  &Locp, &Locq, Atmp, &lld );
            }
            else {
                scalapack_dgebr2d(&icontxt,&scope,&top,
                                  &Locp,&Locq, Atmp, &lld,
                                  &myprow, &jc_proc );
            };
            
            descAtmp[CSRC_] = mypcol;
            /*
             * copy buffer to GPU
             */
            
            inc1 = 1;
            inc2 = 1;
            {
            cublasStatus cu_status = cublasSetVector(isizeAtmp, elemSize,
                                        Atmp, inc1, dAtmp, inc2 );
            CHKERR( cu_status );
            }
            
            // generate Btmp... could be optimized..
            iib = 1;
            jjb = 1;
            descBtmp[RSRC_] = ic_proc;
            if ((transA == 'N') || (transA == 'n')) {
                    transB= 'C';
                    }
            else {
                 transB = 'N';
                 };
            scalapack_pdgeadd( &transB, &k, &n,
                              new_alpha, A, &ia, &ja, descA,
                              new_beta,  Btmp, &iib, &jjb, descBtmp );
            
            Locp = isizeBtmp;
            Locq = 1;
            lld = Locp;
            scope = 'C';
            if (myprow == ic_proc) {
                scalapack_dgebs2d(&icontxt, &scope, &top,
                                  &Locp, &Locq, Btmp, &lld );
            }
            else {
                scalapack_dgebr2d(&icontxt, &scope, &top,
                                  &Locp, &Locq, Btmp, &lld,
                                  &ic_proc, &mypcol );
            };
            
            descBtmp[RSRC_] = myprow;
            /*
             * copy buffer to GPU
             */
            
            inc1 = 1;
            inc2 = 1;
            {
            cublasStatus cu_status = cublasSetVector(isizeBtmp, elemSize,
                                        Btmp, inc1, dBtmp, inc2 );
            CHKERR( cu_status );
            };
        }
        else {
            
            descAtmp[CSRC_] = -1;
            descBtmp[RSRC_] = -1;
            
            iia = 1;
            jja = 1;
            new_alpha = one;
            new_beta = zero;
            

            PROFSTART("syrk:dist_copy_A");
            scalapack_pdgeadd( &transA, &m, &k,
                              new_alpha, A,&ia,&ja,descA,
                              new_beta, Atmp,&iia,&jja,descAtmp );
            PROFEND("syrk:dist_copy_A");
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
            {
            cublasStatus cu_status = cublasSetVector(isizeAtmp, elemSize,
                                        Atmp, inc1, dAtmp, inc2 );
            CHKERR( cu_status );
            }
#endif

       
            
            iib = 1; jjb = 1;
            new_alpha = one;
            new_beta = zero;


            PROFSTART("syrk:dist_copy_B");
            if (use_copy_Atmp) {
              iia = 1;
              jja = 1;
              scalapack_pdgeadd( &transB, &k, &n,
                              new_alpha,  Atmp, &iia, &jja, descAtmp,
                              new_beta, Btmp, &iib, &jjb, descBtmp );
              }
            else {
            if ((transA == 'N') || (transA == 'n')) {
                  transB= 'C';
                  }
            else {
                 transB = 'N';
                 };

            /*
             * B is k by n, it is transpose(A( ia:(ia+n-1), ja:(ja+k-1) ) )
             */
            scalapack_pdgeadd( &transB, &k, &n,
                              new_alpha,  A, &ia, &ja, descA,
                              new_beta, Btmp, &iib, &jjb, descBtmp );

            };
            PROFEND("syrk:dist_copy_B");
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
           {
            cublasStatus cu_status = cublasSetVector(isizeBtmp, elemSize,
                                        Btmp, inc1, dBtmp, inc2 );
            CHKERR( cu_status );
           }
#endif
            
        };
        
        PROFEND("syrk:copy");

#else
/*

//  use GPU direct to copy from A() in GPU device memory 
//  to Dtmp() in  remote GPU device memory on another compute node



        ih = 1; jh = 1;
        scalapack_pdgeadd( notrans, &mm, &nn,
           alpha, A,iia,jja,descA,
           beta,  Dtmp,&ih,&jh,descDtmp );

*/
#endif
        /*
         * Treat as replicated storage
         */
        descAtmp[CSRC_] = mypcol;
        descBtmp[RSRC_] = myprow;

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
         * Do the syrk: 
         * diagonal + upper/lower (depending on the UPLO, 
         *                         upper if 'U', Lower if 'L')
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
            
            zalpha = (double) alpha[REAL_PART];
            zbeta = (double) beta[REAL_PART];



            /*
             * compute large off-diagonal block computed by GEMM
             */
            if ((use_recursion || use_simple) && (m-n >= 1)) {

               if (is_lower) {
                   /*
                    * the triangular diagonal block is at the top
                    * the off-diagonal block is at the bottom
                    */
                   local_extent( m-n,n, ic+n, jc, descC,
                          &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );

                   local_extent( m-n,k, 1+n, 1, descAtmp,
                          &LocpA, &LocqA, &lrA1, &lcA1, &lrA2, &lcA2 );

                   local_extent( k, n, 1, 1, descBtmp,
                          &LocpB, &LocqB, &lrB1, &lcB1, &lrB2, &lcB2 );

                   }
               else {
                   /*
                    * the off-diagonal block is on the top
                    * the triangular diagonal block is at the bottom
                    */
                   local_extent( m-n,n, ic, jc, descC,
                          &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );


                   local_extent( m-n, k, 1,1,descAtmp,
                          &LocpA, &LocqA, &lrA1, &lcA1,  &lrA2, &lcA2 );

                   local_extent( k, n, 1, 1, descBtmp,
                          &LocpB, &LocqB, &lrB1, &lcB1,  &lrB2, &lcB2 );

                   };



                 mm = LocpC;
                 nn = LocqC;
                 kk = k;
                 if ( (mm >= 1) && (nn >= 1) && (kk >= 1) ) {
                      CUBLAS_DGEMM(
                                 CUBLAS_OP_N,CUBLAS_OP_N, mm,nn,kk,
                                 zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                                          (double *) dB(lrB1,lcB1), ldBtmp,
                                 zbeta,   (double *) dC(lrC1,lcC1), ldC );
                      };

               }; /* off-diagonal update */




           if (use_simple) {

              /*
               * loop over column blocks and check for diagonal block to use syrk
               *
               * this may not give high performance since there is a lot of index
               * calculations and the cublas operation are on long slender matrix
               * with width of block size
               *
               * Note take advantage of block cyclic distribution and
               * skip testing some column blocks
               */
              jfinal = jc + n-1;
              jcfirst = Ciafirst( jc, nb, mypcol, descC[CSRC_], npcol );

              for( jstart = jcfirst; jstart <= jfinal; jstart = jend + (npcol-1)*nb + 1) {
                   

                   /* 
                    * make jend at end of block boundary
                    */
                   jend = jstart + (nb - MOD(jstart,nb));
                   if (MOD(jstart,nb) == 0) { jend = jstart; }
                   assert( (jstart <= jend) && (MOD(jend,nb) == 0) );


                   jend = MIN(jend, jfinal);
                   jsize = jend - jstart + 1;

                   assert( mypcol == Cindxg2p( jstart,nb, mypcol, descC[CSRC_], npcol ));

                   /*
                    * check for diagonal block
                    */
                   
                   if (is_lower) {
                     ic_diag = ic + (jstart - jc);
                     jc_diag = jstart;
                     }
                   else {
                     ic_diag = ic + (m-n) + (jstart-jc);
                     jc_diag = jstart;
                     };

                   scalapack_infog2l( &ic_diag, &jc_diag, descC,
                            &nprow,&npcol,  &myprow,&mypcol,
                            &lrC1, &lcC1,   &rsrc,&csrc );

                   /*
                    * this processor owns the diagonal block
                    */
                   has_diagonal = (rsrc == myprow) && (csrc == mypcol);
                   if (has_diagonal) {

                      ia_diag = 1 + (ic_diag - ic );
                      ja_diag = 1;
                      local_extent( jsize,k, ia_diag, ja_diag, descAtmp,
                             &LocpA, &LocqA, &lrA1, &lcA1, &lrA2, &lcA2 );


                      nn = jsize;
                      kk = k;

#ifdef USE_CUBLASV2
                      CUBLAS_DSYRK(
                              (uplo == 'L') || (uplo == 'l') ? 
                                        CUBLAS_FILL_MODE_LOWER :
                                        CUBLAS_FILL_MODE_UPPER,
                             CUBLAS_OP_N,
                             nn, kk, 
                             zalpha, (double *) dA(lrA1,lcA1), ldAtmp,
                             zbeta,  (double *) dC(lrC1,lcC1), ldC );
#else
                      CUBLAS_DSYRK(
                             uplo,  'N', 
                             nn, kk, 
                             zalpha, (double *) dA(lrA1,lcA1), ldAtmp,
                             zbeta,  (double *) dC(lrC1,lcC1), ldC );
#endif



                      }; /* end if (has_diagonal) */


                 /*
                  * process the rest of off-diagonal entries in this column
                  * using GEMM
                  */

                  if (is_lower) {

                     mm = n - (jend - jc+1);
                     nn = jsize;
                     kk = k;

                     iia = 1 + (jend+1-jc);
                     jja = 1;

                     iib = 1;
                     jjb = 1 + (jstart-jc);

                     iic = ic + (jend+1 - jc);
                     jjc = jend;
                     }
                  else {

                     mm = (ic_diag-1) - (ic + (m-n)) + 1;
                     nn = jsize;
                     kk = k;

                     iia = ic + (m-n);
                     jja = 1;

                     iib = 1;
                     jjb = 1 + (jstart-jc);


                     };

                     local_extent( mm, kk, iia, jja, descAtmp,
                         &LocpA, &LocqA, &lrA1, &lcA1, &lrA2, &lcA2 );

                     local_extent( kk, nn, iib, jjb, descBtmp,
                         &LocpB, &LocqB, &lrB1, &lcB1, &lrB2, &lcB2 );

                     local_extent( mm, nn,  iic, jjc, descC,
                         &LocpC, &LocqC, &lrC1, &lcC1, &lrC2, &lcC2 );


                     mm = LocpC;
                     nn = LocqC;
                     kk = k;

                     if ( (mm >= 1) && (nn >= 1) && (kk >= 1)) {

                        assert( LocpC == LocpA );
                        assert( LocqC == LocqB );
                 
#ifdef USE_CUBLASV2
                          CUBLAS_DGEMM( 
                                  CUBLAS_OP_N, CUBLAS_OP_N,
                                   mm,nn,kk,
                                   zalpha,   (double *) dA(lrA1,lcA1), ldAtmp,
                                             (double *) dB(lrB1,lcB1), ldBtmp,
                                   zbeta,    (double *) dC(lrC1,lcC1), ldC );
#else
                          CUBLAS_DGEMM( 'N', 'N', mm,nn,kk,
                                   zalpha,   (double *) dA(lrA1,lcA1), ldAtmp,
                                             (double *) dB(lrB1,lcB1), ldBtmp,
                                   zbeta,    (double *) dC(lrC1,lcC1), ldC );
#endif
                         };



                 }; /* end for(jstart) */

               }
            else if (use_recursion) {

               PROFSTART("Cpsyrk_helper");
               Cpdsyrk_helper( uplo,  n,  k, 
                               alpha,           dAtmp, 1,  1,  descAtmp,
                                                dBtmp, 1,  1,  descBtmp,
                               beta,            C,     ic, jc, descC,  
                               myprow,mypcol,   nprow,npcol );
               PROFEND("Cpsyrk_helper");
                   

               }
       };
#if (0)
            else {
    cublasStatus  cu_status;
    cublasHandle cu_handle; // using the new version of cublas_Dsyrk
    cublasFillMode_t CUDAUPLO;

            
            
            int upper , st_row , st_col, st , en, delta, len_row, dgemm_row;
            upper = (uplo =='U') || (uplo == 'u');
            st = jc + nb * ( ( mypcol - jc_proc + npcol ) % npcol  ); en = jc + m - 1 ;
            st_row = ic + st - jc ;
            
            delta = nb * npcol;
            cublasCreate(&handle);
            
            if ((uplo=='U')||(uplo=='u'))
                CUDAUPLO = CUBLAS_FILL_MODE_UPPER;
            else
                CUDAUPLO = CUBLAS_FILL_MODE_LOWER;
            
            for ( st_col = st ; st_col < en ; st_col = st_col + delta ){
                // Do the syrk on GPU
                local_extent( mb,nb, st_row,st_col,descC,
                             &LocpC, &LocqC, &lrC1,&lcC1, &lrC2,&lcC2 );
                
                jja = 1 ;
                local_extent( mb,k, st_row,jja,descAtmp,
                             &LocpA, &LocqA,
                             &lrA1,&lcA1, &lrA2,&lcA2 );
                iib = 1 ; 
                local_extent( k,nb, iib,st_col,descBtmp,
                             &LocpB, &LocqB,
                             &lrB1,&lcB1, &lrB2,&lcB2 );
                
                /* subroutine of syrk if has work */
                has_work = (LocpC>0)&&(LocqC>0)&&(LocpA>0)&&(LocqA>0)&&(LocpB>0)&&(LocqB>0);
                
                if (has_work) {
                  
                    /* Question
                      1) which subroutine I should use here
                      2) how to make these subroutines run concurrently
                     */
                    // Need to change into col_majored, and others.....!!!!
                    nn = LocpC;
                    kk = k;
                    cublasDsyrk(handle,
                                CUDAUPLO,CUBLAS_OP_N,
                                nn,kk,
                                zalpha,(double *) dA(lrA1,lcA1), ldAtmp,
                                zbeta,(double *) dC(lrC1,lcC1), ldC
                                );
                }
                    
                
                    
                // Do the dgemm on GPU, according to upper or lower
                dgemm_row = st_row + mb ; len_row = ic + m - dgemm_row ;
                if (upper) {len_row = st_row - ic; dgemm_row = ic; }
                
                
                local_extent( len_row,nb, dgemm_row,st_col,descC,
                             &LocpC, &LocqC, &lrC1,&lcC1, &lrC2,&lcC2 );
                
                jja = 1 ;
                local_extent( len_row,k, dgemm_row,jja,descAtmp,
                             &LocpA, &LocqA,
                             &lrA1,&lcA1, &lrA2,&lcA2 );
                iib = 1 ;
                local_extent( k,nb, iib,st_col,descBtmp,
                             &LocpB, &LocqB,
                             &lrB1,&lcB1, &lrB2,&lcB2 );
                
                /* subroutine of dgemm!! */
                has_work = (LocpC>0)&&(LocqC>0)&&(LocpA>0)&&(LocqA>0)&&(LocpB>0)&&(LocqB>0);
                
                if (has_work) {
                    // Need to check !!!
                    mm = LocpC;
                    nn = LocqC;;
                    kk = k;
                    // Question: Is the routine below a non-blocking one? How to make it concurrent if not? Change here?
                    CUBLAS_DGEMM(
                                 CUBLAS_OP_N,CUBLAS_OP_N, mm,nn,kk,
                                 zalpha,  (double *) dA(lrA1,lcA1), ldAtmp,
                                 (double *) dB(lrB1,lcB1), ldBtmp,
                                 zbeta,   (double *) dC(lrC1,lcC1), ldC );
                    
                }
                    
                
                // refresh st_col
                st_row += delta;
            }
            
            cublas_Destroy(handle);
           };
#endif

        
        
    }   
    else {
        /*
         * split operation and use recursion
         */
        k1 = k_max;
        k2 = k - k1;
        
        Cpdsyrk_hhd( uplo, transA, m, n, k1,
                    alpha, A, ia, ja, descA,
                    beta , C, ic, jc, descC );
        
        if (notrans) { iia = ia; jja = ja + k1;}
        else { jja = ja ; iia = ia + k1;};
        
        Cpdsyrk_hhd( uplo, transA, m, n, k2,
                    alpha, A, iia, jja, descA,
                    beta , C, ic, jc, descC );
        
    };
    
    PROFEND("Cpdsyrk_hhd");

}

#ifdef __cplusplus
extern "C"
#endif
void pdsyrk_hhd( char * 	uplo, char  *	transA, int *m, int  * 	n,	int *  	k,
                double *  	alpha, double *  	A, int  * 	ia, int *  	ja, int *  	descA,
                double *  	beta, double *  	C, int  * 	ic, int *  	jc, int *  	descC  )
{
    
    Cpdsyrk_hhd(  *uplo, *transA, *m, *n, *k,
                  alpha, A, *ia, *ja, descA,
                  beta,  C, *ic, *jc, descC );
}


#ifdef __cplusplus
extern "C"
#endif
void pdsyrk_hhd_( char *uplo, char *transA, int *m, int *n, int *k,
                  double *alpha, double *A, int *ia, int *ja, int *descA,
                  double *beta,  double *C, int *ic, int *jc, int *descC )
{
  pdsyrk_hhd( uplo, transA, m,n,k,
              alpha,  A, ia, ja, descA,
              beta,   C, ic, jc, descC );
}
