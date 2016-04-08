#include "ooclu.h"

#if (0)
#define Cpsgemm_hhd(transA,transB,m,n,k, \
    alpha, A,ia,ja,descA, B,ib,jb,descB,beta, C,ic,jc,descC)  \
  scalapack_psgemm(&transA, &transB, &m,&n,&k,  \
      alpha, A, &ia,&ja,descA, B, &ib,&jb,descB, beta, C, &ic,&jc, descC )
#endif

#define ipiv(i)  ipiv_[ IDX1F(i) ]
#define ipiv_hA(i) ipiv_hA_[ IDX1F(i) ]

#define gipiv(i)  gipiv_[IDX1F(i)]

#define dA(i,j)  ( ((float *)A) + IDX2F((i),(j),descA[LLD_]))

#ifdef __cplusplus
extern "C" 
#endif
void psgetrf_gpu(int *m_in, int *n_in, 
   float *A, int *ia_in, int *ja_in, int *descA, 
   int *ipiv_, int *info)
{

int m = *m_in;
int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;

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

int elemSize = sizeof( float );
size_t nbytes;

int nnb, jstart,jend,jsize, isize, jb;
int icontxt, isizeAtmp;


int i,j, iia,jja, ldA, ldhA;
int iinfo = 0;
int iAtmp, jAtmp, iha,jha, iib,jjb,iic,jjc;
int ldAtmp, ldBtmp, lmm,lnn;
int lrA1,lcA1, lrA2,lcA2;

int desc_hA_[DLEN_];
int *desc_hA = &(desc_hA_[0]);

int *ipiv_hA_ = 0;
float *hA = 0;
float *Atmp = 0;
float *dAtmp = 0;

int *gipiv_ = 0;
int desc_Atmp_[DLEN_];
int *desc_Atmp = &(desc_Atmp_[0]);
cublasStatus cu_status;

int isok;
int use_delayed_left_interchange = 1;

int is_mine;
int i1,j1,inc1,  i2,j2,inc2;
int desc_ipiv_hA_[DLEN_];
int *desc_ipiv_hA = &(desc_ipiv_hA_[0]);

int desc_ipiv_[DLEN_];
int *desc_ipiv = &(desc_ipiv_[0]);

int desc_gipiv_[DLEN_];
int *desc_gipiv = &(desc_gipiv_[0]);
int mb,nb, Locp, Locq, lld;


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
/*
 * A is a pointer to GPU device memory but conceptually associated
 * with a scalapack distributed matrix 

 * A is array of complex numbers
 */



*info = 0;

zero[REAL_PART] = 0.0;
zero[IMAG_PART] = 0.0;
one[REAL_PART] = 1.0;
one[IMAG_PART] = 0.0;
neg_one[REAL_PART] = -1.0;
neg_one[IMAG_PART] = 0.0;


/*
 * setup copy of distributed matrix on CPU host
 */

hA = 0;
Atmp = 0;

ictxt = descA[CTXT_];
icontxt = ictxt;

Cblacs_gridinfo( ictxt, &nprow, &npcol,  &myprow, &mypcol );
is_root = (myprow == 0) && (mypcol == 0);
if ((idebug >= 1) && (is_root)) {
  printf("pcgetrf_gpu: m %d n %d ia %d ja %d \n",
      m,n,   ia,ja );
};


ia_proc = Cindxg2p( ia, descA[MB_], myprow, descA[RSRC_], nprow);
ja_proc = Cindxg2p( ja, descA[NB_], mypcol, descA[CSRC_], npcol);


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

  /*
   * setup distribute array hA on host
   */

/*
 * Note, optimal block size on GPU might not be
 * optimal block size on CPU, but assume to be
 * the same for simplicity for now
 */

/*
 * should nnb = descA[NB_] * npcol  ?
 */
nnb = descA[NB_];

minmn = MIN(m,n);
for( jstart=1; jstart <= minmn; jstart = jend + 1) {
  jend = MIN( minmn, jstart + nnb - 1);
  jsize = jend - jstart + 1;

  /*
   * setup matrix on host
   */

  /*
  was iia = (ia-1) + 1;
  */
  j = jstart;
  jb = jsize;

  iia = (ia-1) + jstart;
  jja = (ja-1) + jstart;
  mm = m - jstart + 1;
  nn = jsize;

  if (use_setup_desc) {
    setup_desc( mm,nn, iia,jja,descA, &isize, desc_hA );
    }
  else {
    irsrc = Cindxg2p( iia, descA[MB_], myprow, descA[RSRC_], nprow );
    icsrc = Cindxg2p( jja, descA[NB_], mypcol, descA[CSRC_], npcol );
  
    mb = descA[MB_];
    nb = descA[NB_];
    Locp = Cnumroc( mm, mb, 0,0,nprow );
    Locq = Cnumroc( nn, nb, 0,0,npcol );
    lld = MAX(1,Locp);
    isize = MAX(1,Locp) * MAX(1, Locq );
  
    ictxt = descA[CTXT_];
    iinfo = 0;
    Cdescinit( desc_hA, mm,nn,  mb,nb,  irsrc,icsrc, ictxt, lld, &iinfo);
    assert( iinfo == 0);
    };


  nbytes = isize;
  nbytes *= elemSize;
  if (hA != 0) { 
    free(hA); hA = 0;
  };
  hA = (float *) malloc( nbytes );
  assert( hA != 0 );

  /*
   * distribution of pivot vector is tied to distribution of matrix
   */
  Locp = Cnumroc( desc_hA[M_], desc_hA[MB_], myprow, desc_hA[RSRC_], nprow);
  lld = Locp + desc_hA[MB_];
  nbytes = lld;
  nbytes *= sizeof(int);
  if (ipiv_hA_ != 0) {
    free( ipiv_hA_ ); ipiv_hA_ = 0;
  };
  ipiv_hA_ = (int *) malloc( nbytes );
  assert( ipiv_hA_ != 0);

  Cdescset( desc_ipiv_hA, desc_hA[M_],  1,
              desc_hA[MB_], 1,
              desc_hA[RSRC_], icsrc,
              desc_hA[CTXT_], 
              lld );





  /*
   copy column panel back to CPU host
   to be factored using scalapack
   */ 

  jb = jsize;
  j = jstart;
  mm = m  - j + 1;
  nn = jb;



  /*
    hA(1:mm,1:nn) <-  dA(j:(j+mm-1), j:(j+nn-1) )
   */


  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = 1;
  jj = 1;

  PROFSTART("gpu:hA <- dA");
  Cpsgecopy_d2h( mm,nn, A,iia,jja,descA,  hA, ii,jj, desc_hA );
  PROFEND("gpu:hA <- dA");



  /*
   * factor on host CPU using ScaLAPACK
   * Note the pivot vector is tied to the distribution of the matrix
   * Therefore, we need a different "ipiv_hA" pivot vector
   * that is tied the the distributed matrix hA
   */

  ii = 1;
  jj = 1;
  iinfo = 0;
  mm_lu = mm;
  nn_lu = nn;
  ia_lu = ii;
  ja_lu = jj;

  PROFSTART("gpu:psgetrf");
  scalapack_psgetrf( &mm_lu, &nn_lu, 
        hA, &ia_lu, &ja_lu,  desc_hA, &(ipiv_hA(1)), &iinfo );
  PROFEND("gpu:psgetrf");

  /*
   * broadcast pivot vector to global vector
   */



  i1 = 1;
  j1 = 1;
  inc1 = 1;

  i2 = jstart;
  j2 = 1;
  inc2 = 1;
  mtmp = MIN(mm,nn);
  desc_ipiv_hA[CSRC_] = icsrc;

  use_replicated_storage = FALSE;
  if (use_replicated_storage) {
    int ja_lu_proc;

    ja_lu_proc =   Cindxg2p(ja_lu,desc_hA[NB_],
        mypcol,desc_hA[CSRC_],npcol);

    desc_ipiv_hA[CSRC_] =  ja_lu_proc;

    desc_gipiv[RSRC_] = -1;
    desc_gipiv[CSRC_] = -1;
    scalapack_picopy( &mtmp, &(ipiv_hA(1)), &i1,&j1, desc_ipiv_hA, &inc1,
                        &(gipiv(1)), &i2,&j2, desc_gipiv, &inc2 );
    }
  else {
    /*
     * copy to 1 processors (rsrc,csrc), then
     * broadcast to all processors
     */
        int icontxt = desc_ipiv_hA[CTXT_];
        char scope = 'A'; 
        char top = ' ';
        int ntmp = 1;
        int lld; 

        int ia_lu_proc,ja_lu_proc;
        int rsrc, csrc;

        ia_lu_proc = Cindxg2p( ia_lu, desc_hA[MB_],
               myprow,desc_hA[RSRC_],nprow);
        ja_lu_proc = Cindxg2p( ja_lu, desc_hA[NB_],
               mypcol,desc_hA[CSRC_],npcol);

        rsrc = ia_lu_proc;
        csrc = ja_lu_proc;

        desc_gipiv[RSRC_] = rsrc;
        desc_gipiv[CSRC_] = csrc;
        desc_ipiv_hA[CSRC_] = csrc;

        mtmp = MIN( mm_lu, nn_lu);
        scalapack_picopy( &mtmp, &(ipiv_hA(1)), &i1,&j1,desc_ipiv_hA,&inc1,
                  &(gipiv(1)), &i2,&j2, desc_gipiv, &inc2 );

    if ((myprow == rsrc) && (mypcol == csrc)) {

        lld = mtmp;
        ntmp = 1;
        scalapack_igebs2d( &icontxt, &scope, &top,
            &mtmp, &ntmp, &(gipiv(i2)), &lld );
        }
    else {
      lld = mtmp;
      ntmp = 1;
      scalapack_igebr2d( &icontxt, &scope, &top,
            &mtmp, &ntmp, &(gipiv(i2)), &lld, 
            &rsrc,&csrc );
    };
  };

  if (idebug >= 1) {
    int desctmp[DLEN_];
    char name_ipiv_hA[] = "ipiv_hA";
    char name_gipiv[] = "gipiv";

    if (is_root) {
    printf("jstart %d jend %d \n", jstart,jend);
    printf("mm_lu %d nn_lu %d ia_lu %d ja_lu %d\n",
            mm_lu,   nn_lu,   ia_lu,   ja_lu );
    };

    Cdescset(desctmp, desc_hA[M_], npcol,
        desc_hA[MB_],1,
        desc_hA[RSRC_], desc_hA[CSRC_],
        desc_hA[CTXT_], desc_hA[LLD_] );

    Cpilaprnt( MIN(mm_lu,nn_lu), npcol, &(ipiv_hA(1)), 1,1,desctmp, name_ipiv_hA);

    Cdescset(desctmp, minmn*nprow, npcol,
        minmn, 1,    0,0,
        descA[CTXT_], minmn );
    Cpilaprnt( nprow*minmn, npcol, &(gipiv(1)),1,1,desctmp, name_gipiv);
  };


  /*
   * adjust pivot sequence from 1:min(mm,nn) in ipiv to 
   * jstart:(jstart+min(mm,nn)-1)
   */
    for(int i=1; i <= MIN(mm,nn); i++) {
      i2 = (jstart-1) + i;
      gipiv(i2) = gipiv(i2) + (jstart-1);
    };


  if (iinfo < 0) {
     *info = iinfo;
     return;
     };

  if ((*info == 0) && (iinfo > 0)) {
      *info = iinfo + (j-1);
      return;
      };


  /*
   * transfer factored panel back to GPU device
   */

  iia = (ia-1) + j;
  jja = (ja-1) + j;
  ii = 1;
  jj = 1;
  PROFSTART("gpu:A <- hA");
  Cpsgecopy_h2d(mm,nn, hA, ii,jj, desc_hA,
                       A, iia,jja, descA );
  PROFEND("gpu:A <- hA");





  if (use_delayed_left_interchange) {
    /*
     * do nothing for now
     */
    }
  else {
    /* 
     * apply interchanges to columns 1:(j-1)
     */

    nn = j-1;
    k1 = j;
    k2 = j + jb-1;
    incx = 1;


    PROFSTART("gpu:left swap");
    if (nn >= 1) {
         iia = (ia-1) + 1;
         jja = (ja-1) + 1;
         for(kk=k1; kk <= k2; kk++) {
           ip = gipiv(  kk);
           assert(ip >= kk );
           assert( ip <= m );

           if (kk != ip) {
               inc1 = descA[M_];
               inc2 = descA[M_];
               i1 = (iia-1) + kk;
               i2 = (iia-1) + ip;
               j1 = jja;
               j2 = jja;
               Cpsswap_gpu(nn, A,i1,j1,descA,inc1,
                               A,i2,j2,descA,inc2 );
                };
         };
      };
    PROFEND("gpu:left swap");
    };




  /*
   * apply interchanges to columns (j+jb):n
   */

   nn = n - (jend + 1) + 1;
   k1 = j;
   k2 = j + jb - 1;
   incx = 1;



   PROFSTART("gpu:right swap");
   if (nn >= 1) {
      iia = (ia-1) + 1;
      jja = (ja-1) + (jend+1);
      for(kk=k1; kk <= k2; kk++) {
        ip = gipiv(  kk );
        assert( ip >= kk );
        assert( ip <= m );

        if (ip != kk) {
           i1 = (iia-1) + kk;
           i2 = (iia-1) + ip;
           j1 = jja;
           j2 = jja;
           inc1 = descA[M_];
           inc2 = descA[M_];
           Cpsswap_gpu( nn, A, i1,j1, descA, inc1,
                            A, i2,j2, descA, inc2 );
        };
      };
   };
   PROFEND("gpu:right swap");


   PROFSTART("gpu:pTRSM");


   mm = jb;
   nn = n - (jend+1) + 1;
   if ( (1 <= mm) && (1 <= nn)) {
               /*
               cublasCtrsm('L','L','N','U', mm,nn,
                  alpha, dA(j,j), lddA, dA(j,j+jb), lddA );
               */

     if (use_broadcast_triangular_matrix) {
       /*
        * broadcast triangular part, then solve locally
        */
         char lscope = 'A';
         char ltop = ' ';
         int  msize, nsize, lr1,lc1, lr2,lc2;
         int ia_lu_proc, ja_lu_proc;

       /*
        * copy on local processor
        */

         ia_lu_proc = Cindxg2p(ia_lu, desc_hA[MB_], myprow,
                         desc_hA[RSRC_], nprow );
         ja_lu_proc = Cindxg2p(ja_lu, desc_hA[NB_], mypcol,
                         desc_hA[CSRC_], npcol );

       /*
        * complete mm by mm block on Atmp
        */
       ldAtmp = MAX(1,mm);
       Cdescset(desc_Atmp, mm,mm, mm,mm, 
           ia_lu_proc,ja_lu_proc, icontxt, ldAtmp);
       isizeAtmp = ldAtmp * MAX(1,mm);
       nbytes = isizeAtmp;
       nbytes *= elemSize;

       if (Atmp != 0) { free(Atmp); Atmp = 0; };
       Atmp = (float *) malloc( nbytes );
       assert( Atmp != 0);

#ifdef USE_CUBLASV2
       {
         cudaError_t ierr;
         size_t isize = isizeAtmp;
         isize *= elemSize;

         ierr = cudaMalloc( (void **) &dAtmp, isize );
         assert(ierr == cudaSuccess );
       }
#else
       cu_status = cublasAlloc(isizeAtmp, elemSize, (void **) &dAtmp );
       CHKERR(cu_status);
       assert( dAtmp != 0);
#endif

       ii = 1;
       jj = 1;
       scalapack_psgeadd( notrans, &mm, &mm, 
           one,   hA, &ia_lu, &ja_lu, desc_hA,
           zero,  Atmp, &ii, &jj, desc_Atmp );
                 
       rsrc = desc_Atmp[RSRC_];
       csrc = desc_Atmp[CSRC_];
       if ((myprow == rsrc) && (mypcol == csrc)) {
          scalapack_cgebs2d( &icontxt, &lscope, &ltop,   
              &mm, &mm,  Atmp, &ldAtmp );
          }
       else {
         scalapack_cgebr2d( &icontxt, &lscope, &ltop,
              &mm, &mm, Atmp, &ldAtmp,   &rsrc, &csrc );
       };

       inc1 = 1;
       inc2 = 1;
       cu_status = cublasSetVector(isizeAtmp, elemSize, Atmp, inc1, dAtmp, inc2 );
       CHKERR(cu_status);

       /*
        * perform local solve on GPU
        */
       iia = (ia-1) + j;
       jja = (ja-1) + (j+jb);
       local_extent( mm,nn, iia,jja,descA,  
                    &msize,&nsize, &lr1,&lc1, &lr2,&lc2 );
       if (msize >= 1) {
         assert( msize == mm );
       };

       if ((msize >= 1) && (nsize >= 1)) {
         char lside = 'L';
         char luplo = 'L';
         char ltrans = 'N';
         char ldiag = 'U';

         float zalpha;


         zalpha = (float)1.0;//make_float(1.0,0.0);

         CUBLAS_STRSM( 
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
              mm, nsize, zalpha,
              (float *) dAtmp, ldAtmp,
              dA(lr1,lc1), descA[LLD_] );

       };



       if (Atmp != 0) {
         free(Atmp); Atmp = 0;
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
       CHKERR(cu_status );
#endif


     }
     else {
         /*
          * perform triangular solve using scalapack
          */
         iia = (ia-1) + j;
         jja = (ja-1) + (j+jb);
        setup_desc(mm,nn,iia,jja,descA,  &isize, desc_Atmp );

        nbytes = elemSize;
        nbytes *= isize;
        if (Atmp != 0) {
          free(Atmp); Atmp = 0;
        };
        Atmp = (float *) malloc( nbytes );
        assert( Atmp != 0 );



         /*
          * copy to Atmp(1:mm,1:nn) <- dA(j:(j+mm-1),(j+jb):((j+jb)+nn-1))
          */


         ii = 1; jj = 1;
         PROFSTART("gpu:Atmp <- dA");
         Cpsgecopy_d2h( mm,nn,A,iia,jja,descA,
                           Atmp, ii,jj, desc_Atmp );
         PROFEND("gpu:Atmp <- dA");



         /*
          * perform triangular solve using scalapack
          */

          side = left;
          uplo = lower;
          trans = notrans;
          diag = unit;

          alpha = one;

          iha = 1; 
          jha = 1;
          ii = 1; 
          jj = 1;

          PROFSTART("gpu:pstrsm")
          scalapack_pstrsm( side, uplo, trans, diag, 
              &mm,&nn, alpha,    
              hA, &iha,&jha, desc_hA,
              Atmp,&ii,&jj,  desc_Atmp );
          PROFEND("gpu:pstrsm")
          

          /*
           * copy back to GPU
           */

          iia = (ia-1) + j;
          jja = (ja-1) + (j+jb);
          ii = 1; 
          jj = 1;

          PROFSTART("gpu:A <- Atmp");
          Cpsgecopy_h2d( mm,nn, Atmp,ii,jj,desc_Atmp,
                             A, iia,jja, descA );
          PROFEND("gpu:A <- Atmp");
     };
                           



     };
   PROFEND("gpu:pTRSM");


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
	 cublasSgemm('N','N',mm,nn,kk,
            alpha, dA(j+jb,j),lddA, dA(j,j+jb),lddA,
            beta, dA(j+jb,j+jb), lddA );
         */

        if (use_broadcast_triangular_matrix) {
          /*
           * Copy from GPU to Atmp
           */
          iia = (ia-1) + j;
          jja = (ja-1) + (j+jb);

          setup_desc( kk,nn, iia,jja, descA, &isizeAtmp, desc_Atmp);
          nbytes = isizeAtmp;
          nbytes *= elemSize;
          if (Atmp != 0) { free(Atmp); Atmp = 0; };
          Atmp = (float *) malloc( nbytes );
          assert( Atmp != 0);

          PROFSTART("gpu:Atmp <- A");
          Cpsgecopy_d2h( kk,nn, A,iia,jja,descA, 
                                Atmp,1,1,desc_Atmp );
          PROFEND("gpu:Atmp <- A");
        };


        iic = (ia-1) + (jend+1);
        jjc = (ja-1) + (jend+1);


       iha = jsize+1;
       jha = 1;
       iAtmp = 1; 
       jAtmp = 1;
     

          {
          char transA = 'N';
          char transB = 'N';

          PROFSTART("zgetrf_gpu:psgemm");
          Cpsgemm_hhd( transA, transB, mm,nn,kk, 
           alpha, hA, iha,jha, desc_hA, 
                  Atmp, iAtmp,jAtmp, desc_Atmp, 
           beta,  A, iic,jjc, descA );

          PROFEND("zgetrf_gpu:psgemm");
           };
       };



    if (Atmp != 0) {
       free(Atmp); Atmp = 0;
       };

    if (ipiv_hA_ != 0) {
       free( ipiv_hA_ ); ipiv_hA_ = 0;
       };
    if (hA != 0) {
      free(hA); hA = 0;
      };

   }; /* for (jstart) */


   if (use_delayed_left_interchange) {

     PROFSTART("gpu:dleft swap");
    for(j=1; j <= minmn; j = jend + 1) {
        jend = MIN( minmn, j+nnb-1);
        jsize = jend - j + 1;
        jb = jsize;
        /*
         * apply interchanges to columns 1:(j-1)
         */
   
        nn = j-1;
        k1 = j;
        k2 = j+jb-1;
        incx = 1;
   
   
        if (nn >= 1) {
         iia = (ia-1) + 1; 
         jja = (ja-1) + 1;
         for(kk=k1; kk <= k2; kk++) {
             ip = gipiv(kk);
             assert( ip >= kk );

             if (ip != kk) {
               inc1 = descA[M_];
               inc2 = descA[M_];
               i1 = (iia-1) + kk;
               i2 = (iia-1) + ip;
               j1 = jja;
               j2 = jja;
               Cpsswap_gpu(nn, A, i1,j1,descA, inc1, 
                               A, i2,j2,descA, inc2 );
             };
         };
        };
     }; /* end for j */
     PROFEND("gpu:dleft swap");
   }; /* end if use delayed left interchange */


   /*
    * adjust global pivot from 1:MIN(m,n) to ia:(ia + MIN(m,n)-1)
    * copy global vector back to distributed pivot vector
    */

   for(int j=1; j <= minmn; j++) {
     gipiv(j) = (ia-1) + gipiv(j);
   };


   lld = descA[MB_] + 
         Cnumroc( descA[M_], descA[MB_], myprow, descA[RSRC_], nprow);

   Cdescset( desc_ipiv, 
              descA[M_],1, 
              descA[MB_], 1, 
              descA[RSRC_], -1, descA[CTXT_], lld );

   i1 = 1; j1 = 1; inc1 = 1;
   i2 = ia; j2 = 1; inc2 = 1;
   mtmp = MIN(m,n);

   PROFSTART("gpu:ipiv");
   use_replicated_storage = FALSE;
   if (use_replicated_storage) {
     int msize,nsize,lr1,lc1,lr2,lc2, lrindx,iia;

     local_extent(MIN(m,n),n,ia,ja,descA, &msize,&nsize, &lr1,&lc1, &lr2,&lc2);
     if (msize >= 1) {
       for(lrindx=lr1; lrindx <= lr2; lrindx++) {
         iia = Cindxl2g( lrindx, descA[MB_], myprow, descA[RSRC_], nprow);
         ipiv(lrindx) =  gipiv( (iia-ia) + 1 );
         };
       };
     }
   else  {
     /*
      * copy to a column, then broadcast
      */
     char scope = 'R';
     char top = ' ';
     int Locp, Locq;
     int lld;
     int icontxt = desc_ipiv[CTXT_];

     desc_ipiv[CSRC_] = ja_proc;
     desc_gipiv[RSRC_] = ia_proc;
     desc_gipiv[CSRC_] = ja_proc;

     mtmp = MIN(m,n);
     scalapack_picopy( &mtmp, &(gipiv(1)), &i1,&j1, desc_gipiv, &inc1,
             &(ipiv(1)), &i2, &j2, desc_ipiv, &inc2 );

     if (idebug >= 1) {
       char cmatnm[] = "ipiv after picopy";
       if (is_root) {
         printf("ia_proc %d ja_proc %d i2 %d j2 %d \n",ia_proc,ja_proc,i2,j2);
       };
       Cpilaprnt( mtmp,1, &(ipiv(1)), i2,j2,desc_ipiv, cmatnm);
     };


     Locp = Cnumroc( ia + MIN(m,n)-1, desc_ipiv[MB_], 
                     myprow, desc_ipiv[RSRC_], nprow);
     lld = MAX(1,Locp);
     Locq = 1;
     if (npcol > 1) {
      if (mypcol == ja_proc) {

       scalapack_igebs2d( &icontxt, &scope, &top, 
           &Locp, &Locq,  &(ipiv(1)), &lld );
      }
      else {
       rsrc = myprow;
       scalapack_igebr2d( &icontxt, &scope, &top,
           &Locp, &Locq, &(ipiv(1)), &lld, &rsrc, &ja_proc );
      };
     };

   };
   PROFEND("gpu:ipiv");

     if (idebug >= 1) {
       int desctmp[DLEN_];
       char cmatnm[] = "final ipiv";
       Cdescset( desctmp, 
           descA[M_],npcol,
           descA[MB_],1,
           descA[RSRC_], descA[CSRC_],
           descA[CTXT_], descA[LLD_]);
       Cpilaprnt( MIN(m,n),npcol, &(ipiv(1)), ia,1,desctmp, cmatnm);
     };





  /*
   * clean up
   */
  if (Atmp != 0) {
       free(Atmp); Atmp = 0;
       };
  if (hA != 0) {
       free(hA); hA = 0;
       };
  if (ipiv_hA_ != 0) {
      free( ipiv_hA_ ); ipiv_hA_ = 0;
      };

  if (gipiv_ != 0) {
     free(gipiv_); gipiv_ = 0;
     };


  return;
}



#ifdef __cplusplus
extern "C"
#endif
void psgetrf_gpu_( int *m_in, int *n_in,
           float *A, int *ia_in, int *ja_in, int *descA,
           int *ipiv_,  int *info )
{
    psgetrf_gpu( m_in, n_in,
             (float *) A, ia_in, ja_in, descA,
                ipiv_, info );
}

