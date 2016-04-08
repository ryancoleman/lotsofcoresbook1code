#include "ooclu.h"

#if (0)
#define Cpsgemm_hhd(transA,transB,m,n,k, \
    alpha, A,ia,ja,descA, B,ib,jb,descB,beta, C,ic,jc,descC)  \
  scalapack_psgemm(&transA, &transB, &m,&n,&k,  \
      alpha, A, &ia,&ja,descA, B, &ib,&jb,descB, beta, C, &ic,&jc, descC )
#endif


#define dA(i,j)  ( ((float *)A) + IDX2F((i),(j),descA[LLD_]))

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

void pspotrf_gpu2(char *uplo_in, int *m_in, int *n_in, 
   float *A, int *ia_in, int *ja_in, int *descA, 
   float *Ah, int *iah_in, int *jah_in, int *desc_Ah, 
   int *info)
{

/*
 * uplo_in is ignore for now
 * this routine currently works only on uplo  = "Lower"
 */

float *hA = Ah;


int m = *m_in;
int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;

int iah0 = *iah_in;
int jah0 = *jah_in;


int iproc, jproc;

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

int minmn;
int k1,k2,incx,ip;
int mm, nn, kk, ii, jj, mtmp;
int mm_lu,nn_lu,ia_lu,ja_lu;

int elemSize = sizeof( float );
size_t nbytes;

int nnb, jstart,jend,jsize, isize, jb;
int icontxt, isizeAtmp;


int i,j, iia,jja, ldA, ldhA;
int iinfo = 0;
int has_work = 0;
int iAtmp, jAtmp, iah,jah, iib,jjb,iic,jjc;
int ldAtmp, ldBtmp, lmm,lnn;
int lrA1,lcA1, lrA2,lcA2;




cublasStatus cu_status;

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
char *uplo = lower;
char *trans = notrans;
char *diag = unit;

float zero_[REAL_PART+IMAG_PART+1];
float *zero = &(zero_[0]);

float one_[REAL_PART+IMAG_PART+1];
float *one = &(one_[0]);

float neg_one_[REAL_PART+IMAG_PART+1];
float *neg_one = &(neg_one_[0]);

float beta_[REAL_PART+IMAG_PART+1];
float *beta = &(beta_[0]);

float alpha_[REAL_PART+IMAG_PART+1];
float *alpha = &(alpha_[0]);


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
  printf("pspotrf_gpu2: m %d n %d ia %d ja %d \n",
      m,n,   ia,ja );
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

minmn = MIN(m,n);
for( jstart=1; jstart <= minmn; jstart = jend + 1) {
  jend = MIN( minmn, jstart + nnb - 1);
  jsize = jend - jstart + 1;


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
  Cpsgecopy_d2h_async( mm,nn, A,iia,jja,descA,  Ah, ii,jj, desc_Ah );
  // cublas_sync_stream();
  PROFEND("gpu:hA <- dA");



  /*
   * factor on host CPU using ScaLAPACK
   * (1) factor diagonal block "L11 * L11' = A11"
   * (2) triangular solve   "L21 * L11' = A21"
   *     to get L21
   */
  {

  jb = jsize;
  j = jstart;
  nn = jb;

  /*
   * step (1) Cholesky factor of diagonal block
   */
  

  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = iah0 + (iia-ia);
  jj = jah0 + (jja-ja);





  uplo = lower;
  iinfo = 0;

  PROFSTART("gpu2:pspotrf");
  scalapack_pspotrf( uplo, &nn, Ah, &ii, &jj, desc_Ah, &iinfo );
  PROFEND("gpu2:pspotrf");

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
     PROFSTART("gpu2:pstrsm");
     scalapack_pstrsm( side, uplo, trans, diag, 
               &mm, &nn,  alpha, 
               Ah, &ii, &jj, desc_Ah,
               Ah, &iah, &jah, desc_Ah );
     PROFEND("gpu2:pstrsm");
    };

  }


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
	   cublasSgemm('N','N',mm,nn,kk,
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
            PROFSTART("pspotrf_gpu2:psgemm");

           
           trans = ctranspose;

           scalapack_psgeadd( trans, &kk, &mm, 
              one, 
              hA, &iah, &jah, desc_Ah,
              zero,
              hA, &ii, &jj, desc_Ah );


           /* shiquan debug 09192012 */
           /*if (idebug >= 1) {
             char hA_name3[] = "after transpose hA_1";
             Cpslaprnt( mm,kk,hA,iah,jah,desc_Ah,hA_name3);
             char hA_name4[] = "after transpose hA_2";
             Cpslaprnt( kk,mm,hA,ii,jj,desc_Ah,hA_name4);  
           } */


            Cpsgemm_hhd( transA, transB, mm,nn,kk, 
               alpha, hA, iah,jah, desc_Ah, 
                  hA, ii,jj,   desc_Ah, 
               beta,  A, iic,jjc, descA );
              PROFEND("pspotrf_gpu2:psgemm");
           }
          else {
            char c_uplo = 'L';
            char c_trans = 'N';

              PROFSTART("pspotrf_gpu2:pssyrk");
            Cpssyrk_hhd( c_uplo, c_trans, mm,nn,kk, 
            alpha, hA, iah,jah, desc_Ah, 
            beta,  A, iic,jjc, descA );
              PROFEND("pspotrf_gpu2:pssyrk");
           };
      };
   }; /* end for (jstart) */

  return;
}



#ifdef __cplusplus
extern "C"
#endif
void pspotrf_gpu2_( char *uplo, int *m_in, int *n_in,
           float *A, int *ia_in, int *ja_in, int *descA,
           float *Ah, int *iah_in, int *jah_in, int *desc_Ah,
           int *info )
{
    pspotrf_gpu2( uplo, m_in, n_in,
             (float *) A, ia_in, ja_in, descA,
             (float *) Ah, iah_in, jah_in, desc_Ah,
                 info );
}

