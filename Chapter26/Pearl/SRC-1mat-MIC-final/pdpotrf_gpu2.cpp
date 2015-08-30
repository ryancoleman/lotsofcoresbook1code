#include "ooclu.h"


#if (0)
#define Cpdgemm_hhd(transA,transB,m,n,k, \
    alpha, A,ia,ja,descA, B,ib,jb,descB,beta, C,ic,jc,descC)  \
  scalapack_pdgemm(&transA, &transB, &m,&n,&k,  \
      alpha, A, &ia,&ja,descA, B, &ib,&jb,descB, beta, C, &ic,&jc, descC )
#endif


#define dA(i,j)  ( ((double *)A) + IDX2F((i),(j),descA[LLD_]))

#ifdef __cplusplus
extern "C" 
#endif
/*
 * Compute Cholesky factorization for matrix A(ia,ja) on GPU
 * Note this is a rectangular m by n panel
 *
 * This is the image of Ah(iah,jah) on CPU
 * The final factorization will be on CPU as well
 */

void pdpotrf_gpu2(char *uplo_in, int *m_in, int *n_in, 
   double *A, int *ia_in, int *ja_in, int *descA, 
   double *Ah, int *iah_in, int *jah_in, int *desc_Ah, 
   int *info)
{

/*
 * uplo_in is ignore for now
 * this routine currently works only on uplo  = "Lower"
 */

double *hA = Ah;


int m = *m_in;
int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;

int iah0 = *iah_in;
int jah0 = *jah_in;


int iproc, jproc;


#ifndef USE_MIC
cublasStatus cu_status;
#endif

/*
 * set use_gemm = FALSE to use the more efficient
 * SYRK to exploit symmetry
 */
const int use_gemm = FALSE;

const int idebug = 1;

int ia_proc, ja_proc;
int lrindx, lcindx, rsrc,csrc, irsrc,icsrc;
int ictxt, nprow,npcol, myprow,mypcol;
int is_root;

int lr1,lc1,lr2,lc2;

int minmn;
int k1,k2,incx,ip;
int mm, nn, kk, ii, jj, mtmp;
int mm_lu,nn_lu,ia_lu,ja_lu;

int elemSize = sizeof( double );
size_t nbytes;

int nnb, jstart,jend,jsize, isize, jb;
int icontxt, isizeAtmp;


int i,j, iia,jja, ldA, ldhA;
int iinfo = 0;
int has_work = 0;
int iAtmp, jAtmp, iah,jah, iib,jjb,iic,jjc;
int ldAtmp, ldBtmp, lmm,lnn;
int lrA1,lcA1, lrA2,lcA2;





int isok;
int use_delayed_left_interchange = 1;

int is_mine;
int i1,j1,inc1,  i2,j2,inc2;


int mb,nb, Locp, Locq, lld;



char direc = 'F';
char rowcol = 'R';

char upper[] = "Upper";
char lower[] = "Lower";

char left[] = "Left";
char right[] = "Right";

char transpose[] = "Transpose";
char ctranspose[] = "ConjuateTranspose";
char notrans[] = "NoTrans";

char unit[] = "Unit";
char nonunit[] = "NonUnit";

char *side = right;
char *uplo = uplo_in;
char *trans = notrans;
char *diag = unit;

int is_lower = (*uplo == 'L') || (*uplo == 'l');

double zero_[REAL_PART+IMAG_PART+1];
double *zero = &(zero_[0]);

double one_[REAL_PART+IMAG_PART+1];
double *one = &(one_[0]);

double neg_one_[REAL_PART+IMAG_PART+1];
double *neg_one = &(neg_one_[0]);

double beta_[REAL_PART+IMAG_PART+1];
double *beta = &(beta_[0]);

double alpha_[REAL_PART+IMAG_PART+1];
double *alpha = &(alpha_[0]);


const int use_gpu_factor = FALSE;
int descDtmp[DLEN_];
double *Dtmp = 0;
double *hDtmp = 0;
int ih, jh;

/* test code */
double diagonal_block_trace1, diagonal_block_trace2;
/* test code */

/*
 * A is a pointer to GPU device memory but conceptually associated
 * with a scalapack distributed matrix 

 * A is array of complex numbers in GPU device memory
 * Ah is array in CPU memory
 */



*info = 0;

zero[REAL_PART] = 0.0;
zero[IMAG_PART] = 0.0;
one[REAL_PART] = 1.0;
one[IMAG_PART] = 0.0;
neg_one[REAL_PART] = -1.0;
neg_one[IMAG_PART] = 0.0;




ictxt = descA[CTXT_];
icontxt = ictxt;

Cblacs_gridinfo( ictxt, &nprow, &npcol,  &myprow, &mypcol );
is_root = (myprow == 0) && (mypcol == 0);
if ((idebug >= 1) && (is_root)) {
  printf("pdpotrf_gpu2: m %d n %d mn %f ia %d ja %d \n",
      m,n, float(m)*float(n), ia,ja );
};


ia_proc = Cindxg2p( ia, descA[MB_], myprow, descA[RSRC_], nprow);
ja_proc = Cindxg2p( ja, descA[NB_], mypcol, descA[CSRC_], npcol);


/*
 * Note, optimal block size on GPU might not be
 * optimal block size on CPU, but assume to be
 * the same for simplicity for now
 */

/*
 * should nnb = descA[NB_] * npcol  ?
 * or nnb = descA[NB_]
 */

nnb = descA[NB_];


/*
 * outer Main loop 
 */

if (use_gpu_factor) {
 /*
  * setup data structures
  */

  int mbnb = MAX( descA[NB_], descA[MB_] );
  size_t nbytes = 0;
  int isizeDtmp = 0;
  int iinfo = 0;
  


  isizeDtmp = mbnb*mbnb;
  nbytes = isizeDtmp;
  nbytes *= elemSize;
  #ifdef USE_MIC
    hDtmp = (double*) malloc(nbytes);
    assert( hDtmp != 0);
    Dtmp = (double*) offload_Alloc(nbytes);
    assert( Dtmp != 0);
  #else
    cudaError_t cuda_err = cudaMallocHost( &hDtmp, nbytes );
    assert( cuda_err == cudaSuccess );
 
    CUBLAS_MALLOC( Dtmp,  isizeDtmp, elemSize );
  #endif

  Cdescinit(descDtmp, mbnb, mbnb,
           mbnb,mbnb, myprow,mypcol, descA[CTXT_],mbnb,&iinfo);
  assert( iinfo == 0);

};
  
  


minmn = MIN(m,n);
for( jstart=1; jstart <= minmn; jstart = jend + 1) {
  int j,jb;

  jend = MIN( minmn, jstart + nnb - 1);
  jsize = jend - jstart + 1;
  j = jstart;
  jb = jsize;
  

  if (use_gpu_factor) {


/* test code, compare between cpu llt decompose diagonal block and gpu one */
  /*copy column panel back to CPU host to be factored using scalapack */ 
  /*Ah(ii:(ii+mm-1),jj:(jj+nn-1)) <-  dA(j:(j+mm-1), j:(j+nn-1) ) */
  jb = jsize;
  j = jstart;
  mm = m  - j + 1;
  nn = jb;
  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = iah0 + (iia-ia);
  jj = jah0 + (jja-ja);
  PROFSTART("gpu:hA <- dA");
  Cpdgecopy_d2h_async( mm,nn, A,iia,jja,descA,  Ah, ii,jj, desc_Ah );
  // cublas_sync_stream();
  PROFEND("gpu:hA <- dA");
   /* factor on host CPU using ScaLAPACK
   * (1) factor diagonal block "L11 * L11' = A11" */
  uplo = lower;
  iinfo = 0;
  PROFSTART("gpu2:pdpotrf");
  scalapack_pdpotrf( uplo, &nn, Ah, &ii, &jj, desc_Ah, &iinfo );
  PROFEND("gpu2:pdpotrf");
  if (iinfo < 0) {
    *info = iinfo;
     return;
    };
  if ((*info == 0) && (iinfo > 0)) {
     *info = iinfo + (j-1);
     };
  diagonal_block_trace1 = 0.0;
  for( i=1; i <= nn; i++) { 
    diagonal_block_trace1 += Ah[ii+i,jj+i];
  }
/* test code, compare between cpu llt decompose diagonal block and gpu one */



   /*
    *  (1) perform Cholesky factorization on GPU,
    *  (2) broadcast factor to processors in same processor column
    *  (3) perform triangular solve
    *  (4) copy processed column to CPU
    */

   /*
    * (1) perform Cholesky factorizatioin on GPU
    */


    iia = (ia-1) + jstart;
    jja = (ja-1) + jstart;

    int iia_proc = Cindxg2p( iia, descA[MB_], myprow, descA[RSRC_],nprow);
    int jja_proc = Cindxg2p( jja, descA[NB_], mypcol, descA[CSRC_],npcol);

    is_mine = (iia_proc == myprow) && (jja_proc == mypcol);

    if (is_mine) {
      int ierr = 0;

      local_extent( jsize,jsize, iia,jja,descA, 
         &mm, &nn, 
         &lr1, &lc1,  &lr2, &lc2 );


       assert( mm == nn );
       assert( nn >= 1 );
       assert( mm == jsize );




       iinfo = 0;
       PROFSTART("gpu2:potrf_gpu");
       #ifdef USE_MIC
         offload_dpotrf(uplo, &nn, dA(lr1,lc1), &descA[LLD_], &iinfo);
       #else
       ierr = magma_dpotrf_gpu( *uplo, nn, dA(lr1,lc1), 
                 descA[LLD_],  &iinfo);
       #endif
       PROFEND("gpu2:potrf_gpu");

       if (iinfo != 0) {
          printf("dpotrf return iinfo=%d \n", iinfo);
          printf("nn %d lr1 %d lc1 %d \n", nn,lr1,lc1 );
          };



      if (iinfo < 0) {
        *info = iinfo;
         return;
        };

      if ((*info == 0) && (iinfo > 0)) {
         *info = iinfo + (j-1);
         return;
         };
       };

      /*
       * (2) copy and broadcast  factor to all processors in
       * same processor column
       */
      mm = jsize;
      nn = jsize;
      iia = (ia-1) + j;
      jja = (ja-1) + j;
      ii = iah0 + (iia-ia);
      jj = jah0 + (jja-ja);
      Cpdgecopy_d2h( mm,nn, A,iia,jja,descA,  Ah, ii,jj, desc_Ah );

      /*
       * set replicated storage to perform broadcast
       */
      descDtmp[CSRC_] = jja_proc;
      descDtmp[RSRC_] = -1;

      PROFSTART("gpu2:bcast hDtmp");
      mm = jsize;
      nn = jsize;
      alpha = one;
      beta = zero;
      ih = 1;
      jh = 1;
      scalapack_pdgeadd( notrans, &mm, &nn, 
           alpha, Ah, &ii,&jj,desc_Ah,
           beta,  hDtmp,&ih,&jh,descDtmp );
      PROFEND("gpu2:bcast hDtmp");

      /*
       * copy to GPU
       */
      descDtmp[RSRC_] = myprow;
      descDtmp[CSRC_] = mypcol;

      Cpdgecopy_h2d( mm,nn,
             hDtmp,1,1,descDtmp,
             Dtmp,1,1,descDtmp );

      /*
      -- debug
      printf("after h2d: hDtmp %lf \n",
            *hDtmp  );
      printf("is_lower %d *uplo=%c \n",is_lower, *uplo );
      */
      
      /*
       * perform triangular solve
       * solve L21 * L11' = A21
       */
      if (is_lower) {
       int LocpA, LocqA;

       mm = m - (jend+1) + 1;
       nn = jsize;
       iia = (ia-1) + (jend+1);
       jja = (ja-1) + jstart;
       local_extent( mm,nn, iia,jja,descA, 
         &LocpA, &LocqA, 
         &lr1, &lc1,  &lr2, &lc2 );

       /*
       -- debug
       printf("debug: m %d n %d jstart %d  jend %d mm %d nn %d\n",
                      m,   n,   jstart,    jend,   mm,   nn );
       printf("iia %d jja %d \n", iia,jja);
       printf("debug: LocpA %d LocqA %d lr1 %d lc1 %d\n",
              LocpA,LocqA,lr1,lc1 );
       */

       has_work = (LocpA >= 1) && (LocqA >= 1);
       if (has_work) {
          /*
           solve L21 * L11 = A21 
           */

         char lside = 'R';
         char luplo = 'L';
         char ltrans = 'C';
         char ldiag = 'N';

           mm = LocpA;
           nn = LocqA;

           alpha = one;


         PROFSTART("gpu2:trsm");
         #ifdef USE_MIC
                 offload_dtrsm(
                      &lside, &luplo, &ltrans, &ldiag, 
                      &mm, &nn, 
                      alpha,
                      Dtmp, &descDtmp[LLD_],
                      dA(lr1,lc1), &descA[LLD_] );
         #else
           #ifdef USE_CUBLASV2
           { 
                /*
                 CUBLAS_DTRSM(
                     ((lside == 'l')||(lside == 'L')) ?
                        CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT,
                     ((luplo == 'l')||(luplo == 'L')) ?
                        CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER,
                     ((ltrans == 'c')||(ltrans == 'C')) ?
                       CUBLAS_OP_C :
                         ((ltrans == 't')||(ltrans == 'T')) ?
                            CUBLAS_OP_T : CUBLAS_OP_N,
                     ((ldiag == 'u')||(ldiag == 'U')) ?
                        CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT,
                      mm, nn, 
                       alpha,
                       Dtmp, descDtmp[LLD_],
                      dA(lr1,lc1), descA[LLD_] );
                  */
                magmablas_dtrsm( lside, luplo, ltrans, ldiag,
                    mm,nn,
                    *alpha, Dtmp, descDtmp[LLD_],
                    dA(lr1,lc1), descA[LLD_] );
              }
            #else
                 CUBLAS_DTRSM(
                      lside, luplo, ltrans, ldiag, 
                      mm, nn, 
                      *alpha,
                      Dtmp, descDtmp[LLD_],
                      dA(lr1,lc1), descA[LLD_] );
            #endif
          #endif
          PROFEND("gpu2:trsm");

          };

         

         
        }
      else {
        /*
         * upper triangular not implemented yet
         */
       };
       
       /*
        -- debug
       printf("debug: after step (3) \n");
       */

       /*
        * (4) copy processed column to CPU
        */

       if (is_lower) {
          
        mm = m - (jend+1) + 1;
        nn = jsize;

        iia = (ia-1) + (jend+1);
        jja = (ja-1) + (jstart);
        ii = iah0 + (iia - ia);
        jj = jah0 + (jja - ja);
        alpha = one;
        beta = zero;

        Cpdgecopy_d2h( mm, nn, 
           A,  iia, jja, descA,
           hA, ii,  jj,  desc_Ah );

          }
       else {
         /*
          * upper triangular not implemented yet
          */
        };

       /*
        -- debug
       printf("debug: after step (4) \n");
       */


/* test code, compare between cpu llt decompose diagonal block and gpu one */
  diagonal_block_trace2 = 0.0;
  for( i=1; i <= nn; i++) { 
    diagonal_block_trace2 += Ah[ii+i,jj+i];
  }
  if  (  pow( (diagonal_block_trace2-diagonal_block_trace1),2.0 )  >1.0E-8) {
    printf("pdpotrf_gpu2: diagonal_block_trace1, diagonal_block_trace2=%f, %f\n", diagonal_block_trace1, diagonal_block_trace2);
  }
/* test code, compare between cpu llt decompose diagonal block and gpu one */

  }
  else {
  /*
   copy column panel back to CPU host
   to be factored using scalapack
   */ 



  /*
    Ah(ii:(ii+mm-1),jj:(jj+nn-1)) <-  dA(j:(j+mm-1), j:(j+nn-1) )
   */


  jb = jsize;
  j = jstart;
  mm = m  - j + 1;
  nn = jb;

  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = iah0 + (iia-ia);
  jj = jah0 + (jja-ja);



  PROFSTART("gpu:hA <- dA");
  Cpdgecopy_d2h_async( mm,nn, A,iia,jja,descA,  Ah, ii,jj, desc_Ah );
  // cublas_sync_stream();
  PROFEND("gpu:hA <- dA");



  /*
   * factor on host CPU using ScaLAPACK
   * (1) factor diagonal block "L11 * L11' = A11"
   * (2) triangular solve   "L21 * L11' = A21"
   *     to get L21
   */

  /*
   * step (1) Cholesky factor of diagonal block
   */
  

  uplo = lower;
  iinfo = 0;

  PROFSTART("gpu2:pdpotrf");
  {
//  translate the coordinates and call pdpotrf with i,j = 1
//  fix some problems to
//  scalapack_pdpotrf( uplo, &nn, Ah, &ii, &jj, desc_Ah, &iinfo );
//  when local matrix is large (~40000)
    double *pAh = 0;
    int mAh = nn;
    int nAh = nn;
    int mbAh = desc_Ah[MB_];
    int nbAh = desc_Ah[NB_];
    int descpAh[DLEN_];
    int p_one = 1;
    int rsrcpAh = Cindxg2p( ii, desc_Ah[MB_], myprow, desc_Ah[RSRC_],nprow);
    int csrcpAh = Cindxg2p( jj, desc_Ah[NB_], mypcol, desc_Ah[CSRC_],npcol);
    int pia = Cnumroc( ii-1, mbAh, myprow, desc_Ah[RSRC_], nprow);
    int pja = Cnumroc( jj-1, nbAh, mypcol, desc_Ah[CSRC_], npcol);
    Cdescinit(descpAh, mAh, nAh, mbAh, nbAh,
        rsrcpAh, csrcpAh, ictxt, desc_Ah[LLD_], &iinfo);
    pAh = ((double *)Ah) + IDX2F(pia+1,pja+1,descpAh[LLD_]);
//  printf("gpu2:pdpotrf skip no\n");
    scalapack_pdpotrf( uplo, &nn, pAh, &p_one, &p_one, descpAh, &iinfo );
  } 
//scalapack_pdpotrf( uplo, &nn, Ah, &ii, &jj, desc_Ah, &iinfo );
  PROFEND("gpu2:pdpotrf");

  if (iinfo < 0) {
    *info = iinfo;
     return;
    };

  if ((*info == 0) && (iinfo > 0)) {
     *info = iinfo + (j-1);
     };

  /*
   * step (2) triangular solve "L21 * L11' = A21 "
   */

  side = right;
  uplo = lower; 
  trans = ctranspose;
  diag = nonunit;

  nn = jsize;
  mm = m - (jend + 1) + 1;
  alpha = one;
  iah = ii + nn;
  jah = jj;

  has_work = (mm >= 1) && (nn >= 1);
  if (has_work) {
     PROFSTART("gpu2:pdtrsm");
//   printf("gpu2:pdtrsm skip no\n");
     scalapack_pdtrsm( side, uplo, trans, diag, 
               &mm, &nn,  alpha, 
               Ah, &ii, &jj, desc_Ah,
               Ah, &iah, &jah, desc_Ah );
     PROFEND("gpu2:pdtrsm");
    };

 };


    /* 
     * -------------------------
     * update trailing submatrix
     * -------------------------
     */


	alpha = neg_one;
	beta = one;
	mm = m-(jend+1) + 1;
	nn = n-(jend+1) + 1;
	kk = jb;

 
      has_work = (1 <= mm) && (1 <= nn) && (1 <= kk);
      if (has_work) {
        
        /*
	     cublasDgemm('N','N',mm,nn,kk,
            alpha, dA(j+jb,j),lddA, dA(j,j+jb),lddA,
            beta, dA(j+jb,j+jb), lddA );
         */


          char transA = 'N';
          char transB = 'N';

          iah = (iah0-1) + (j+jb);
          jah = (jah0-1) + j;

          ii = (iah0-1) + j;
          jj = (jah0-1) + (j+jb);

          iic = (ia-1) + (j+jb);
          jjc = (ja-1) + (j+jb);



          if (use_gemm) {
            PROFSTART("pdpotrf_gpu2:pdgemm");
             PROFSTART("trans_copyA");
           trans = ctranspose;
           scalapack_pdgeadd( trans, &kk, &mm, 
              one, 
              hA, &iah, &jah, desc_Ah,
              zero,
              hA, &ii, &jj, desc_Ah );
             PROFEND("trans_copyA");

            Cpdgemm_hhd( transA, transB, mm,nn,kk, 
               alpha, hA, iah,jah, desc_Ah, 
                  hA, ii,jj,   desc_Ah, 
               beta,  A, iic,jjc, descA );
              PROFEND("pdpotrf_gpu2:pdgemm");
           }
          else {
            char c_uplo = 'L';
            char c_trans = 'N';

              PROFSTART("pdpotrf_gpu2:pdsyrk");
            Cpdsyrk_hhd( c_uplo, c_trans, mm,nn,kk, 
            alpha, hA, iah,jah, desc_Ah, 
            beta,  A, iic,jjc, descA );
              PROFEND("pdpotrf_gpu2:pdsyrk");

           };

      };



   }; /* end for (jstart) */

if (use_gpu_factor) {
  /*
   * clean up
   */
  #ifdef USE_MIC
   free(hDtmp);
   offload_Free(Dtmp);
  #else
   cudaError_t   cuda_err;


   cuda_err = cudaFreeHost( (void *) hDtmp );
   assert( cuda_err == cudaSuccess );

   CUBLAS_FREE( Dtmp );
  #endif

   };

  return;
}



#ifdef __cplusplus
extern "C"
#endif
void pdpotrf_gpu2_( char *uplo, int *m_in, int *n_in,
           double *A, int *ia_in, int *ja_in, int *descA,
           double *Ah, int *iah_in, int *jah_in, int *desc_Ah,
           int *info )
{
    pdpotrf_gpu2( uplo, m_in, n_in,
             (double *) A, ia_in, ja_in, descA,
             (double *) Ah, iah_in, jah_in, desc_Ah,
                 info );
}

