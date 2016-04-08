#include "ooclu.h"

#if (0)
#define Cpcgemm_hhd(transA,transB,m,n,k, \
    alpha, A,ia,ja,descA, B,ib,jb,descB,beta, C,ic,jc,descC)  \
  scalapack_pcgemm(&transA, &transB, &m,&n,&k,  \
      alpha, A, &ia,&ja,descA, B, &ib,&jb,descB, beta, C, &ic,&jc, descC )
#endif

#define ipiv(i)  ipiv_[ IDX1F(i) ]
#define ipiv_Ah(i) ipiv_Ah_[ IDX1F(i) ]
#define lpvt(i) lpvt_[ IDX1F(i) ]


#define gipiv(i)  gipiv_[IDX1F(i)]

#define dA(i,j)  ( ((cuComplex *)A) + IDX2F((i),(j),descA[LLD_]))

#ifdef __cplusplus
extern "C" 
#endif
/*
 * Compute LU factorization for matrix A(ia,ja) on GPU
 * This is the image of Ah(iah,jah) on CPU
 * The final factorization will be on CPU as well
 */

void pcgetrf_gpu2(int *m_in, int *n_in, 
   cuComplex *A, int *ia_in, int *ja_in, int *descA, 
   float *Ah, int *iah_in, int *jah_in, int *desc_Ah, 
   int *ipiv_Ah_, int *info)
{

float *hA = Ah;


int m = *m_in;
int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;

int iah0 = *iah_in;
int jah0 = *jah_in;


int iproc, jproc;



const int use_split_gemm = FALSE;
const int use_delayed_left_swap = TRUE;

const int use_setup_desc = TRUE;
const int idebug = 0;
int use_replicated_storage = FALSE;
const int use_broadcast_triangular_matrix = TRUE;

int ia_proc, ja_proc;
int lrindx, lcindx, rsrc,csrc, irsrc,icsrc;
int ictxt, nprow,npcol, myprow,mypcol;
int is_root;

int minmn;
int k1,k2,incx,ip;
int mm, nn, kk, ii, jj, mtmp;
int mm_lu,nn_lu,ia_lu,ja_lu;

int elemSize = sizeof( cuComplex );
size_t nbytes;

int nnb, jstart,jend,jsize, isize, jb;
int icontxt, isizeAtmp;


int i,j, iia,jja, ldA, ldhA;
int iinfo = 0;
int iAtmp, jAtmp, iah,jah, iib,jjb,iic,jjc;
int ldAtmp, ldBtmp, lmm,lnn;
int lrA1,lcA1, lrA2,lcA2;


int desc_ipiv_Ah_[DLEN_];
int *desc_ipiv_Ah = &(desc_ipiv_Ah_[0]);

int desc_lpvt_[DLEN_];
int *desc_lpvt = &(desc_lpvt_[0]);

cublasStatus cu_status;

int *gipiv_ = 0;
int isok;
int use_delayed_left_interchange = 1;

int is_mine;
int i1,j1,inc1,  i2,j2,inc2;

int desc_ipiv_[DLEN_];
int *desc_ipiv = &(desc_ipiv_[0]);

int desc_gipiv_[DLEN_];
int *desc_gipiv = &(desc_gipiv_[0]);
int mb,nb, Locp, Locq, lld;

int *lpvt_ = 0;


char direc = 'F';
char rowcol = 'R';

char left[] = "Left";
char lower[] = "Lower";
char notrans[] = "NoTrans";
char unit[] = "Unit";

char *side = left;
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


short *has_copied = 0;
/*
 * A is a pointer to GPU device memory but conceptually associated
 * with a scalapack distributed matrix 

 * A is array of complex numbers in GPU device memory
 * Ah is array in CPU memory
 */

has_copied = (short *) malloc( sizeof(short) * m+1 );
assert(has_copied != 0 );
for( int i=0; i <= m; i++) {
  has_copied[i] = FALSE;
};


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
  printf("pcgetrf_gpu2: m %d n %d ia %d ja %d \n",
      m,n,   ia,ja );
};


ia_proc = Cindxg2p( ia, descA[MB_], myprow, descA[RSRC_], nprow);
ja_proc = Cindxg2p( ja, descA[NB_], mypcol, descA[CSRC_], npcol);

#if (0)
/*
 * setup global pivot vector
 */
lld = MIN(m,n) + descA[MB_];
nbytes = lld;
nbytes *= sizeof(int);
if (gipiv_ != 0) {
  free(gipiv_); gipiv_ = 0;
};
gipiv_ = (int *) malloc( nbytes );
assert( gipiv_ != 0 );


desc_gipiv[DTYPE_] = descA[DTYPE_];
desc_gipiv[CTXT_] = descA[CTXT_];
desc_gipiv[M_] = MIN(m,n);
desc_gipiv[N_] = 1;
desc_gipiv[MB_] = desc_gipiv[M_];
desc_gipiv[NB_] = desc_gipiv[N_];
desc_gipiv[LLD_] = lld;

desc_gipiv[RSRC_] = -1;
desc_gipiv[CSRC_] = -1;
#endif

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

lpvt_ = (int *) malloc( sizeof(int) * nnb );
assert( lpvt_ != 0 );

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
  Cpcgecopy_d2h_async( mm,nn, A,iia,jja,descA,  Ah, ii,jj, desc_Ah );
  // cublas_sync_stream();
  PROFEND("gpu:hA <- dA");



  /*
   * start copy of "U" part
   */
  {
  PROFSTART("d2h:U part");
   j = jstart;
   jb = jsize;
   int mm = jsize;
   int nn = n - (jend + 1) + 1;

   int iia = (ia-1) + j;
   int jja = (ja-1) + j+jb;
   int iah = iah0 + (iia-ia);
   int jah = jah0 + (jja-ja);

   Cpcgecopy_d2h_async( mm,nn,A,iia,jja,descA,
         hA, iah, jah, desc_Ah );
  PROFEND("d2h:U part");
  }


  /*
   * factor on host CPU using ScaLAPACK
   * Note the pivot vector is tied to the distribution of the matrix
   * Therefore, we need a different "ipiv_Ah" pivot vector
   * that is tied the the distributed matrix hA
   */
  {
  jb = jsize;
  j = jstart;
  mm = m  - j + 1;
  nn = jb;

  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = iah0 + (iia-ia);
  jj = jah0 + (jja-ja);

  iinfo = 0;
  mm_lu = mm;
  nn_lu = nn;
  ia_lu = ii;
  ja_lu = jj;

  PROFSTART("gpu:pcgetrf");
  scalapack_pcgetrf( &mm_lu, &nn_lu, 
        Ah, &ia_lu, &ja_lu,  desc_Ah, &(ipiv_Ah(1)), &iinfo );
  PROFEND("gpu:pcgetrf");



  if (iinfo < 0) {
     *info = iinfo;
     return;
     };

  if ((*info == 0) && (iinfo > 0)) {
      *info = iinfo + (j-1);
      return;
      };
  }





  /*
   * apply interchanges to columns (j+jb):n
   * by
   * copying  A( j:jend, (j+jb):n ) to CPU
   * copying  A( ipvt(k1:k2), (j+jb):n) to CPU 
   * perform  swap on CPU, then copy  A( ipvt(k1:k2), (j+jb):n) 
   * back to GPU
   */

   PROFSTART("gpu:right swap");

   j = jstart;
   jb = jsize;
   mm = jsize;
   nn = n - (jend + 1) + 1;

   iia = (ia-1) + j;
   jja = (ja-1) + j+jb;
   iah = iah0 + (iia-ia);
   jah = jah0 + (jja-ja);

   /*
   Cpcgecopy_d2h_async( mm,nn,A,iia,jja,descA,
         hA, iah, jah, desc_Ah );
    */


   /*
    * broadcast  part of pivot vector
    * use the last MB entries of pivot vector
    */
   iah = (iah0-1) + j;
   jah = (jah0-1) + j;
   scalapack_infog2l( &iah, &jah, desc_Ah, &nprow, &npcol, &myprow, &mypcol, 
         &lrindx, &lcindx, &iproc, &jproc );


   Locp = Cnumroc( desc_Ah[M_], desc_Ah[MB_], myprow, desc_Ah[RSRC_], nprow );

   if (jsize > desc_Ah[MB_]) {

     int inc1 = 1;
     int inc2 = 1;
     int i1 = (iah0-1) + jstart;
     int j1 = 1;
     int i2 = 1;
     int j2 = 1;

     Cdescset(  desc_ipiv_Ah, 
         desc_Ah[M_] + desc_Ah[MB_],1,     desc_Ah[MB_], 1, 
         desc_Ah[RSRC_], desc_Ah[CSRC_], 
         desc_Ah[CTXT_], Locp + desc_Ah[MB_]);
     desc_ipiv_Ah[CSRC_] = jproc;


     /*
      * collect to processor (iproc,jproc)
      */
     Cdescset( desc_lpvt, nnb,1,  nnb,1,  iproc,jproc, 
               desc_Ah[CTXT_], nnb );

     scalapack_picopy( &jsize, 
           &(ipiv_Ah(1)), &i1, &j1, desc_ipiv_Ah, &inc1,
           &(lpvt(1)), &i2, &j2, desc_lpvt, &inc2 );

     /*
      * broadcast to all
      */
       if ((myprow == iproc) && (mypcol == jproc)) {
         char scope[] = "All";
         char top[] = " ";
         int m0 = jsize;
         int n0 = 1;
         int lld = m0;

         scalapack_igebs2d( &icontxt, scope, top, &m0, &n0, 
           &(lpvt(1)),  &lld);
         }
       else {
         char scope[] = "All";
         char top[] = " ";
         int m0 = jsize;
         int n0 = 1;
         int lld = m0;

         scalapack_igebr2d( &icontxt, scope, top, &m0, &n0, 
            &(lpvt(1)), &lld, &iproc, &jproc);
       };


      }
     else
   {

     // char scope[] = "All";
     char scope[] = "Col";
     char top[] = " ";

   // if ((myprow == iproc)  && (mypcol == jproc)) {
   if ((myprow == iproc)  ) {
     int m0, n0, lld; 

     for(i=1; i <= jsize; i++) {
       lpvt( i ) = ipiv_Ah( (lrindx-1) + i );
       };

     m0 = jsize;
     n0 = 1;
     lld = m0;

     /*
      * use broadcast "All" for simplify, can consider using only "Column"
      */

       scalapack_igebs2d( &icontxt, scope, top, &m0, &n0, 
         &(lpvt(1)),  &lld);

   }
   else {

     int m0,n0,lld;

     m0 = jsize;
     n0 = 1;
     lld = m0;


       scalapack_igebr2d( &icontxt, scope, top, &m0, &n0, 
         &(lpvt(1)), &lld, &iproc, &jproc);
     };

   };

   /*
    * determine which rows to copy
    */

   PROFSTART("pcgetrf_gpu2:d2h for swap");
   {
     int irowmax = (iah0-1) + jend;
     for (i=1; i <= jsize; i++) {
      int ip = lpvt(i); 
      int is_outside_block = (ip > irowmax);
      int irow = ip - iah0;

      if (is_outside_block) {
        iia =  ia + (ip-iah0);
        jja = (ja-1) + (j+jb);

        iah = ip;
        jah = (jah0-1) + (j+jb);
        if (!has_copied[irow]) {
          Cpcgecopy_d2h( 1, nn, A, iia,jja,descA, 
               hA, iah, jah, desc_Ah );
          has_copied[irow] = TRUE;
          };
       };
     };
     /*
      * reset boolean vector
      */
     for( int i=1; i <= jsize; i++) {
       int ip = lpvt(i); 
       int irow = ip  - iah0;
       has_copied[irow] = FALSE;
       };
          

   }
   PROFEND("pcgetrf_gpu2:d2h for swap");

   /* 
    * perform row swap on CPU
    */

   // cublas_sync_stream();
   PROFSTART("pcgetrf_gpu2: swap on CPU");
   {
     int inc1 = desc_Ah[M_]; 
     int inc2 = desc_Ah[M_];
     int jah_right = (jah0-1) + (jend+1);
     int jah_left = (jah0-1) + 1;
     int nn_right = n - (jend+1)+1;
     int nn_left = (j-1);

     for(int i=1; i <= jsize; i++) {

       int ip = lpvt(i); 
       int irow = (iah0-1) + (j-1) + i;

       assert( (1 <= irow) && (irow <= desc_Ah[M_]));
       assert( (1 <= ip) && (ip <= desc_Ah[M_]) );
       assert( ip >= irow);

       if (ip != irow) {
         if (nn_right >= 1) {
           scalapack_pcswap(&nn_right, 
               Ah, &irow, &jah_right, desc_Ah, &inc1,
               Ah, &ip,   &jah_right, desc_Ah, &inc2);
         };
         if ((!use_delayed_left_swap) && (nn_left >= 1) ) {
           scalapack_pcswap(&nn_left,  
               Ah, &irow, &jah_left, desc_Ah, &inc1,
               Ah, &ip,   &jah_left, desc_Ah, &inc2);
         };
       };
     };
   }
   PROFEND("pcgetrf_gpu2: swap on CPU");



   /*
    * copy rows back to GPU
    */



   PROFSTART("pcgetrf_gpu2:h2d for swap");
   {
     int irowmax =  (iah0-1) + jend;
     for (i=1; i <= jsize; i++) {
      int ip = lpvt(i); 
      int is_outside_block = (ip > irowmax);
      int irow = ip - iah0;

      if (is_outside_block) {
       iia =  ia + (ip-iah0);
       jja = (ja-1) + (j+jb);

       iah = ip;
       jah = (jah0-1) + (j+jb);
       if (!has_copied[irow]) {
         Cpcgecopy_h2d_async( 1, nn,  hA,iah,jah, desc_Ah, A, iia,jja,descA);
         has_copied[irow] = TRUE;
         };
       };
     };

     /* 
      * reset boolean vector
      */
     for(int i=1; i <= jsize; i++) {
       int ip = lpvt(i); 
       int irow = ip - iah0;
       has_copied[irow] = FALSE;
     };

     // cublas_sync_stream();

   };
   PROFEND("pcgetrf_gpu2:h2d for swap");



   PROFEND("gpu:right swap");







         /*
          * perform triangular solve using scalapack
          */
       {

          side = left;
          uplo = lower;
          trans = notrans;
          diag = unit;

          alpha = one;

          j = jstart;
          iah = (iah0-1) + j; 
          jah = (jah0-1) + j;
          ii = iah; 
          jj = (jah0-1) + (jend+1);

          mm = jsize;
          nn = n - (jend+1) + 1;

          PROFSTART("gpu:pctrsm")
          scalapack_pctrsm( side, uplo, trans, diag, 
              &mm,&nn, alpha,    
              hA, &iah,&jah, desc_Ah,
              hA, &ii, &jj, desc_Ah );
          PROFEND("gpu:pctrsm")
         }



    /*
     * update trailing submatrix
     */


	alpha = neg_one;
	beta = one;
	mm = m-(jend+1) + 1;
	nn = n-(jend+1) + 1;
	kk = jb;

 
      if ((1 <= mm) && (1 <= nn) && (1 <= kk)) {
        
        /*
	     cublasCgemm('N','N',mm,nn,kk,
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

          PROFSTART("zgetrf_gpu2:pcgemm");
          if (use_split_gemm) {
            /*
             * split operation to compute next panel
             */
            int mm1 = mm;
            int jstart2 = jend + 1;
            int jend2 = MIN(n, jstart2 + nnb-1);
            int nn1 = jend2 - jstart2 + 1;
            int kk1 = kk;

            if ((mm1 >= 1) && (nn1 >= 1) && (kk1 >= 1)) {
              Cpcgemm_hhd( transA,transB,mm1,nn1,kk1,
                alpha, hA,iah,jah,desc_Ah,
                       hA,ii,jj,  desc_Ah,
                beta,  A,iic,jjc, descA );
              };


            int ii2 = ii;
            int jj2 = jj + nn1;
            int iic2 = iic;
            int jjc2 = jjc + nn1;
            int mm2 = mm;
            int nn2 = nn - nn1;
            int kk2 = kk;

            if ((mm2 >= 1) && (nn2 >= 1) && (kk2 >= 1)) {
              Cpcgemm_hhd( transA,transB,mm2,nn2,kk2,
                alpha, hA,iah,jah,desc_Ah,
                       hA,ii2,jj2,desc_Ah,
                beta,  A,iic2,jjc2,descA );
                };

          }
          else {
          Cpcgemm_hhd( transA, transB, mm,nn,kk, 
           alpha, hA, iah,jah, desc_Ah, 
                  hA, ii,jj,   desc_Ah, 
           beta,  A, iic,jjc, descA );
          };

          PROFEND("zgetrf_gpu2:pcgemm");
           };


      /*
       * perform left swap
       */


   if (use_delayed_left_swap) {
     PROFSTART("zgetrf_gpu2:left swap");

     int inc1 = desc_Ah[M_]; 
     int inc2 = desc_Ah[M_];
     int jah_right = (jah0-1) + (jend+1);
     int jah_left = (jah0-1) + 1;
     int nn_right = n - (jend+1)+1;
     int nn_left = (j-1);

     for(int i=1; i <= jsize; i++) {

       int ip = lpvt(i); 
       int irow = (iah0-1) + (j-1) + i;

       assert( (1 <= irow) && (irow <= desc_Ah[M_]));
       assert( (1 <= ip) && (ip <= desc_Ah[M_]) );
       assert( ip >= irow);

       if (ip != irow) {
         if (nn_left >= 1) {
           scalapack_pcswap(&nn_left,  
               Ah, &irow, &jah_left, desc_Ah, &inc1,
               Ah, &ip,   &jah_left, desc_Ah, &inc2);
         };
       };
     };
   PROFEND("zgetrf_gpu2:left swap");
   };




   }; /* for (jstart) */


  if (has_copied != 0) {
     free( has_copied );
  };
  if (lpvt_ != 0) {
    free(lpvt_);
  };

  return;
}



#ifdef __cplusplus
extern "C"
#endif
void pcgetrf_gpu2_( int *m_in, int *n_in,
           float *A, int *ia_in, int *ja_in, int *descA,
           float *Ah, int *iah_in, int *jah_in, int *desc_Ah,
           int *ipiv_,  int *info )
{
    pcgetrf_gpu2( m_in, n_in,
             (cuComplex *) A, ia_in, ja_in, descA,
             (float *) Ah, iah_in, jah_in, desc_Ah,
                ipiv_, info );
}

