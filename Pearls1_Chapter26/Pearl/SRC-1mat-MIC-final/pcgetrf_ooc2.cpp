#include "ooclu.h"


#define ipiv(i) ipiv_[IDX1F(i)]

#define ipivY(i) ipivY_[IDX1F(i)]
#define gipiv(i) gipiv_[IDX1F(i)]

#ifdef __cplusplus
extern "C" 
#endif
void pcgetrf_ooc2( int *m_in, int *n_in, 
             float *A, int *ia_in, int *ja_in, int *descA, 
	     int *ipiv_,  int *memsize_in, int *info )
{
int m = *m_in;
int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;
int memsize = *memsize_in;

int use_replicated_storage = FALSE;
const int use_pcgetrf_gpu2 = TRUE;

const int idebug = 0;
const int use_immediate_copy_to_host = TRUE;


char left[] = "Left";
char lower[] = "Lower";
char unit[] = "Unit";
char notrans[] = "NoTrans";

char *side = left;
char *uplo = lower;
char *trans = notrans;
char *diag = unit;

char forward[] = "Forward";
char rowwise[] = "Rowwise";
char *direct = forward;
char *rowcol = rowwise;

cublasStatus cu_status;

int nprow,npcol,myprow,mypcol, icontxt;
int is_root;


float *dX = 0;
cuComplex *dY = 0;

int mb,nb, nbx, nby, npanels, isize, isizeX,isizeY;
int ix0,jx0, iix,jjx, ix,iy, jx,jy;
int descY[DLEN_];

int descX_[DLEN_];
int *descX = &(descX_[0]);

int k1,k2,incx,ip;
int jstart,jend,jsize,i,j,jb, jstartm1;
int kfinal, kstart,kend,ksize, knb, ihY,jhY;
int isizehY;
int has_work;

int mm_lu,nn_lu,iy_lu,jy_lu;
int msize, minmn, mm,nn, kk, iia,jja, ii,jj, lld;
int iproc, rsrc, csrc, Locp, Locq, irsrc,icsrc;

int elemSize = sizeof(cuComplex);

float *hY = 0;
int *ipivY_ = 0;
int deschY[DLEN_];
size_t nbytes;
int i1,j1,inc1,  i2,j2,inc2, lrindx;

int desc_ipivY_[DLEN_];
int *desc_ipivY = &(desc_ipivY_[0]);

int desc_ipiv_[DLEN_];
int *desc_ipiv = &(desc_ipiv_[0]);

int desc_gipiv_[DLEN_];
int *desc_gipiv = &(desc_gipiv_[0]);
int *gipiv_;

const int use_new_hY = FALSE;
const int do_adjust_pivot = FALSE;
int is_mine;

int ia_proc,ja_proc;

int isok;
int iinfo = 0;
char direc = 'F';
char row = 'R';

float beta_[REAL_PART+IMAG_PART+1];
float *beta = &(beta_[0]);

float alpha_[REAL_PART+IMAG_PART+1];
float *alpha = &(alpha_[0]);

float one_[REAL_PART+IMAG_PART+1];
float *one = &(one_[0]);

float zero_[REAL_PART+IMAG_PART+1];
float *zero = &(zero_[0]);

float neg_one_[REAL_PART+IMAG_PART+1];
float *neg_one = &(neg_one_[0]);


zero[REAL_PART] = 0.0;
zero[IMAG_PART] = 0.0;

one[REAL_PART] = 1.0;
one[IMAG_PART] = 0.0;

neg_one[REAL_PART] = -1;
neg_one[IMAG_PART] = 0.0;


dX = 0;
dY = 0;
hY = 0;

*info = 0;

/*
 * perform LU factorization similar to scalapack PCGETRF
 * but use memsize entries on GPU
 */



/*
 * check arguments
 */

if ( (m < 0) || (m > descA[M_])) {
   *info = -1;
   return;
   };

if ((n < 0) || (n > descA[N_])) {
  *info = -2;
  return;
  };

if ((ia < 1) || (ia > descA[M_])) {
   *info = -4;
   return;
   };

if ((ja < 1) || (ja > descA[N_])) {
   *info = -5;
   return;
   };

if (memsize <= 0) {
  *info = -8;
  return;
  };


/*
 * estimate storage
 * for panels X and Y
 * width of panel X is nbx
 * width of paney Y is nby
 */

icontxt = descA[CTXT_];
Cblacs_gridinfo( icontxt, &nprow,&npcol, &myprow, &mypcol );
is_root = (myprow == 0) && (mypcol == 0);


/*
 * Note that optimal block size for GPU may not be
 * optimal block size for CPU.
 * Assume they are the same for simplicity for now.
 */

mb = descA[MB_];
nb = descA[NB_];


/*
 * should nbx be larger?
   nbx = descA[NB_]*npcol;
 */
nbx = descA[NB_];

/*
 * estimate storage for a global panel of mm by (nb*npcol)
 */


mm = m;
nn = nb*npcol;

ia_proc = Cindxg2p(ia,mb, myprow, descA[RSRC_], nprow );
ja_proc = Cindxg2p(ja,nb, mypcol, descA[CSRC_], npcol );

rsrc = ia_proc;
csrc = ja_proc;

/*
 * Slightly over estimate
 * storage and to have the same and consistent
 * result across all processors
 */
Locp = Cnumroc( mm, mb, 0, 0, nprow );
Locq = Cnumroc( nn, nb, 0, 0, npcol );
isize = MAX(1, Locp)*MAX(1,Locq );



npanels = MAX(1,memsize/isize);
nby = npanels * nn;
nby = MAX( nbx, nby );

if ((is_root)  && (idebug >= 0)) {
  printf("pcgetrf_ooc2: mm %d nn %d mb %d nb %d isize %d \n",
                        mm,nn, mb,nb, isize );
  printf("pcgetrf_ooc2: memsize %d npanels %d nby %d\n",
      memsize, npanels, nby );
};

minmn = MIN(m,n);


/*
 * directly use host matrix as "X" panel
 * no need to allocate storage for "X" panel
 */

 dX = A;
 ix = ia;
 jx = ja;
 descX = descA;

 /*
  * setup replicated global pivot vector gipiv(1:MIN(m,n))
  * each processor has a complete gipiv vector
  */
 lld = minmn + descA[MB_];
 nbytes = lld;
 nbytes *= sizeof(int);
 gipiv_ = (int *) malloc( nbytes );
 Cdescset( desc_gipiv,
     minmn, 1,
     minmn, 1,
     -1,-1,
     descA[CTXT_], lld );

 /*
  * Preallocate storage for Y panel and pivot vector
  */
 setup_desc( m,nby, ia,ja,descA,  &isizeY, descY );
 nbytes = isizeY;
 nbytes *= elemSize;

#ifdef USE_CUBLASV2
{
  cudaError_t ierr;
  size_t isize = isizeY;
  isize *= elemSize;

  ierr = cudaMalloc( (void **) &dY, isize );
  assert(ierr == cudaSuccess );
}
#else
 cu_status = cublasAlloc( isizeY, elemSize, (void **) &dY );
 CHKERR(cu_status);
#endif

 lld = descY[LLD_] + descY[MB_];
 Cdescset( desc_ipivY,
     descY[M_], 1,
     descY[MB_], 1,
     descY[RSRC_], -1,
     descY[CTXT_], lld );

 nbytes = lld;
 nbytes *= elemSize;
 ipivY_ = (int *) malloc( nbytes );



 /*
  * Main outer loop over panels
  */
 for(jstart=1; jstart <= n; jstart = jend + 1 ) {
    j = jstart;
    jend = MIN( n, jstart + nby - 1 );
    jsize = jend - jstart + 1;
    jb = jsize;

    /*
     copy Y(1:m,1:jsize) <-  A(1:m,jstart:jend)
     */

     iy = 1; 
     jy = 1;
     iia = (ia-1) + 1;
     jja = (ja-1) + jstart;
     mm = m;
     nn = jsize;

     PROFSTART("Y <- A");
     Cpcgecopy_h2d( mm,nn, A,iia,jja,descA,   dY, iy,jy,descY );
     PROFEND("Y <- A");


     kfinal = MIN(m, jstart-1);
     knb = nbx;


     if (use_new_hY) {
       /*
        * use a separate region in host memory
        * otherwise, use the original "A" matrix itself
        */
       setup_desc(nbx,nby, 1,1,descY,  &isizehY, deschY );
       if (hY != 0) { free(hY); hY = 0; };
       nbytes = isizehY;
       nbytes *= elemSize;
       hY = (float *) malloc( nbytes );
       assert( hY != 0 );
     }
     else {
       /*
        * reuse part of original "A" matrix as hY
        * copy descriptor
        */
       hY = A;
       memcpy( &(deschY[0]), &(descA[0]), sizeof(deschY) );
     }

     for(kstart=1; kstart <= kfinal; kstart = kend + 1) {
             kend = MIN(kfinal, kstart+knb-1);
             ksize = kend - kstart + 1;


/*
 *             -----------------------------------------------
 *             d_X(kstart:m,1:ksize) <=> A(kstart:m,kstart:kend)
 *             -----------------------------------------------
 */


              /*
               * directly use host matrix as X panel
               * X(kstart,1) is equivalent to A(iia,jja)
               * iia = (ia-1) + kstart, jja = (ja-1) + jstart
               */

              dX = A;
              descX = descA;

/*
 *            ------------------ 
 *            update part of d_Y
 *
 *            use scalapack to perform triangular solve
 *
 *            d_X(kstart:kend,kstart:kend)\d_Y(kstart:kend,1:jsize)
 *
 *            need to copy dY(kstart:kend,1:jsize) to hY(1:ksize,1:jsize)
 *            ------------------ 
 */
              PROFSTART("X:Ctrsm");

              mm = ksize;
              nn = jsize;
              iy = kstart;
              jy = 1;
              if (use_new_hY) {
                setup_desc(mm,nn, iy,jy,descY,  &isize, deschY );
                assert( isize <= isizehY );


                ihY = 1;
                jhY = 1;
              }
              else {
                /*
                 * reuse part of original "A" matrix
                 */
                hY = A;
                ihY = (ia-1) + iy;
                jhY = (ja-1) + (jstart-1) + jy;
              };
              Cpcgecopy_d2h( mm,nn, dY, iy,jy,   descY, 
                                    hY, ihY,jhY, deschY );
              cublas_sync_stream();

              /*
                 Ctrsm('L','L','N','U', mm,nn,alpha,
                  X(kstart,1), ldx, Y(kstart,1), ldy );
               */

               iix = (ia-1) + kstart;
               jjx = (ja-1) + kstart;
               if (use_new_hY) {
                 ihY = 1;
                 jhY = 1;
                 }
               else {
                 ihY = (ia-1) + kstart;
                 jhY = (ja-1) + jstart;
               };


              if (idebug >= 1) {
                int irprnt = 0;
                int icprnt = 0;
                int nout = 6;
                int ilen;
                char hY_name[] = "hY";
                char dX_name[] = "dX";


                if ((myprow==0) && (mypcol==0)) {
                  printf("=== jstart %d jend %d  kstart %d kend %d \n", 
                              jstart,   jend,    kstart,   kend );
                };

                Cpclaprnt( mm,nn,hY,ihY,jhY,descY,hY_name);

                /*
                   { 
                   float * work;
                   size_t nbytes;
                   nbytes  = sizeof(float)*2*descY[MB_]*2;
                   nbytes = MAX(nbytes, 1024*1024);

                   work = (float *) malloc(  nbytes  );
                   ilen = strlen(hY_name);
                   scalapack_pclaprnt( &mm,&nn,hY,&ihY,&jhY,deschY,
                    &irprnt,&icprnt,
                    hY_name, &nout, work, &ilen );
                    free( work );
                    }
                    */

                Cpclaprnt(mm,mm,dX,iix,jjx,descX,dX_name);

                /*
                   {
                   float *work;
                   size_t nbytes;
                   nbytes =  sizeof(float)*2*descX[MB_]*2;
                   nbytes = MAX(nbytes, 1024*1024 );

                   work = (float *) malloc( nbytes );
                   ilen = strlen(dX_name);
                   scalapack_pclaprnt( &mm,&mm, dX,&iix,&jjx,descX,
                     &irprnt,&icprnt,
                    dX_name, &nout, work, &ilen );
                    free(work);
                    }
                    */
              };



               alpha = one;
               scalapack_pctrsm(side,uplo,trans,diag,
                           &mm,&nn, alpha,
                           dX,&iix,&jjx,descX,
                           hY,&ihY,&jhY,deschY );

              PROFEND("X:Ctrsm");




/*
 *            ------------------------
 *            update lower part of d_Y
 *            k1 = kstart + mm
 *            d_Y(k1:m,1:jsize) <- d_Y(k1:m,1:jsize) -  
 *                d_X(k1:m,ksize) * d_Y(kstart:kend, 1:jsize)
 *            -------------------------
 */
             PROFSTART("X:Cgemm");


             alpha = neg_one;
             beta = one;
             k1 = kend+1;
             mm = m - k1 + 1;
             nn = jsize;
             kk = ksize;

             has_work = (mm >= 1) && (nn >= 1) && (kk >= 1);
             if (has_work)  {
              iix = (ia-1) + k1;
              jjx = (ja-1) + kstart;
              iy = k1;
              jy = 1;
              if (use_new_hY) {
                ihY = 1;
                jhY = 1;
                }
              else {
                ihY = (ia-1) + kstart;
                jhY = (ja-1) + (jstart-1) + 1;
              };

              if (idebug >= 1) {

                char cmatnm[] = "dX1";

                int irprnt = 0;
                int icprnt = 0;
                int nout = 6;
                int ilen;

                ilen = strlen(cmatnm);

                if ((myprow==0) && (mypcol == 0)) {
                  printf("mm %d kk %d iix %d jjx %d \n", 
                          mm,   kk,   iix,   jjx );
                };

                Cpclaprnt(mm,kk,dX,iix,jjx,descX,cmatnm);

                /*
                   {
                   float *work;
                   size_t nbytes;

                   nbytes = sizeof(float)*2*descX[MB_]*2;
                   nbytes = MAX( nbytes, 1024*1024 );
                   work = (float *) malloc( nbytes );
                  scalapack_pclaprnt( &mm,&kk, dX, &iix,&jjx,descX, 
                    &irprnt,&icprnt, cmatnm, &nout, work, &ilen );
                    free(work);
                    }
                    */
              };


              Cpcgemm_hhd('N', 'N', mm,nn,kk,
                          alpha, dX,iix,jjx,descX,
                                 hY,ihY,jhY,deschY,
                          beta,  dY, iy,jy,descY );

              };
              PROFEND("X:Cgemm");



              if (use_immediate_copy_to_host) {
                /*
                 * copy results back to CPU, 
                 * since it is not used on GPU anyway
                 * A(kstart:kend), jstart:jend) <- hA
                 */
                if (use_new_hY) {
                 alpha = one;
                 beta = zero;
                 iia = (ia-1) + kstart;
                 jja = (ja-1) + jstart;
                 ihY = 1;
                 jhY = 1;

                 mm = ksize;
                 nn = jsize;

                 scalapack_pcgeadd( notrans, &mm,&nn,
                     alpha, hY, &ihY,&jhY, deschY, 
                     beta,  A,  &iia,&jja, descA );
                   };
                 
                 }
              else {
              /*
               * copy results back to GPU
               */
                iy = kstart;
                jy = 1;
                if (use_new_hY) {
                  ihY = 1;
                  jhY = 1;
                  }
                else {
                  ihY = (ia-1) + iy;
                  jhY = (ja-1) + (jstart-1) + jy;
                };

                mm = ksize;
                nn = jsize;
                Cpcgecopy_h2d( mm,nn, hY,ihY,jhY,deschY,
                             dY, iy,jy,descY );
                 };




             }; /* end for kstart */

            if (use_new_hY) {
              if (hY != 0) {
                  free(hY); hY = 0;
                 };
            };

/*
 *           -------------------------------------
 *           all left updates applied to panel d_Y
 *           perform LU factorization of  panel Y 
 *           in GPU device memory
 *           -------------------------------------
 */
            j = jstart;
            jsize = jend - jstart + 1;
            jb  = jsize;
            mm = m - j + 1;
            nn = jsize;

             /*
              * similar to pcgetrf( &mm,&nn,Y,j,1,&descY, &(ipivY(1)), &iinfo);
             */
             iy = j; 
             jy = 1;
             has_work = (mm >= 1) && (nn >= 1);
             if (has_work) {
               iinfo = 0;
               mm_lu = mm;
               nn_lu = nn;
               iy_lu = iy;
               jy_lu = jy;


               if (use_pcgetrf_gpu2) {
                 int iah0 = (ia-1) + iy_lu;
                 int jah0 = (ja-1) + jstart;

                 PROFSTART("Y:pcgetrf_gpu2");
                 pcgetrf_gpu2( &mm_lu, &nn_lu, 
                   dY,&iy_lu,&jy_lu,descY, 
                   A, &iah0, &jah0, descA,
                   &(ipiv(1)), &iinfo );
                 PROFEND("Y:pcgetrf_gpu2");

                 /* copy ipiv(:) to ipivY(:) */
                    {
                      int msizeA,nsizeA,lr1A,lc1A,lr2A,lc2A;
                      int msizeY,nsizeY,lr1Y,lc1Y,lr2Y,lc2Y;

                      local_extent( mm_lu,1, iah0,jah0,descA,
                          &msizeA,&nsizeA,
                          &lr1A,&lc1A, &lr2A,&lc2A );
                      local_extent( mm_lu,1, iy_lu,jy_lu,descY,
                          &msizeY,&nsizeY,
                          &lr1Y,&lc1Y, &lr2Y,&lc2Y );

                      assert(msizeA == msizeY);
                      for(int i=1; i <= msizeA; i++) {
                        ipivY( (lr1Y-1) + i ) = ipiv( (lr1A-1) + i );
                        };
                    }
                 }
               else {
                 PROFSTART("Y:pcgetrf_gpu");
                 pcgetrf_gpu( &mm_lu, &nn_lu, 
                   dY,&iy_lu,&jy_lu,descY, &(ipivY(1)), &iinfo );
                 PROFEND("Y:pcgetrf_gpu");
               };
             };

#ifdef USE_FAKE_CUBLAS
             if (idebug >= 2) {
               char cmatnm[] = "dY2";
               if (is_root) { printf("after pcgetrf_gpu\n"); };
               Cpclaprnt( mm,nn, (float *) dY,iy,jy,descY,cmatnm);
             };
#endif

             if (has_work) {
               if (iinfo < 0) {
                 *info = iinfo;
                 if (idebug >= 1) {
                   printf("pcgetrf_ooc2: pcgetrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
               if (iinfo > 0) {
                 *info = iinfo + (j-1);
                 if (idebug >= 1) {
                   printf("pcgetrf_ooc2: pcgetrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
             };


/*
 *            -------------------------------------
 *            A(1:m,jstart:jend) = Y(1:m,1:jsize)
 *            -------------------------------------
 */

            if (!use_pcgetrf_gpu2) {
            PROFSTART("A <- Y");


            if (use_immediate_copy_to_host) {
              /*
               * Just copy the rectangular submatrix
               * that was just factored
               * Y(jstart:m, jstart:jend) 
               */
              j = jstart;
              mm = m - j + 1;
              nn = jsize;

              iy = j; 
              jy = 1;
              iia = (ia-1) + jstart;
              jja = (ja-1) + jstart;
              if ((mm >= 1) && (nn >= 1)) {
               Cpcgecopy_d2h( mm,nn, dY,iy,jy,descY,
                              A,iia,jja,descA );
                };
              }
            else {
              /*
               * transfer the entire Y panel
               *
               * A (ia:(ia+m-1),jstart:jend) <- Y(1:m,1:jsize)
               */

              mm = m;
              nn = jsize;
              iia = (ia-1) + 1;
              jja = (ja-1) + jstart;
              Cpcgecopy_d2h( mm,nn,  dY,1,1,descY, 
                                   A,iia,jja,descA );
            };

            PROFEND("A <- Y");
            };


            /*
             * broadcast pivot vector to gipiv
             */
            i1 = jstart;
            j1 = 1;
            inc1 = 1;
            i2 = jstart;
            j2 = 1;
            inc2 = 1;

            mm = m - jstart + 1;
            nn = jsize;
            msize = MIN(mm,nn);

            PROFSTART("Y:gipiv");

            use_replicated_storage = FALSE;
            if (use_replicated_storage) {

            desc_ipivY[CSRC_] = descY[CSRC_];
            desc_gipiv[RSRC_] = -1;
            desc_gipiv[CSRC_] = -1;


            if (msize >= 1) {
             scalapack_picopy( &msize,  &(ipivY(1)), &i1,&j1, desc_ipivY, &inc1,
                 &(gipiv(1)), &i2,&j2, desc_gipiv, &inc2 );
             };
                  
            }
            else {
              char scope = 'A';
              char top = ' ';
              int lld;

              int iy_lu_proc,jy_lu_proc;

              iy_lu_proc = Cindxg2p( iy_lu, 
                  descY[MB_],myprow,descY[RSRC_],nprow);
              jy_lu_proc = Cindxg2p( jy_lu,
                  descY[NB_],mypcol,descY[CSRC_],npcol);
              
              rsrc = iy_lu_proc;
              csrc = jy_lu_proc;

              desc_ipivY[CSRC_] = csrc;

              desc_gipiv[RSRC_] = rsrc;
              desc_gipiv[CSRC_] = csrc;

              if (msize >= 1) {
               scalapack_picopy( &msize,  
                 &(ipivY(1)), &i1,&j1, desc_ipivY, &inc1,
                 &(gipiv(1)), &i2,&j2, desc_gipiv, &inc2 );
              };

              /*
               * broadcast to All
               */

              scope = 'A';
              top = ' ';
              nn = 1;
              lld = msize;
              if ((myprow == rsrc) && (mypcol == csrc)) {
                if (msize >= 1) {
                  scalapack_igebs2d(&icontxt, &scope, &top, 
                       &msize, &nn, &(gipiv(i2)), &lld );
                  };
              }
              else {
                if (msize >= 1) {
                  scalapack_igebr2d(&icontxt, &scope, &top, 
                       &msize, &nn, &(gipiv(i2)), &lld,   &rsrc,&csrc );
                 };

              };

            };
            PROFEND("Y:gipiv");
/*
 *           ------------------------------------------
 *           apply interchanges to columns 1:(jstart-1)
 *           ------------------------------------------
 */

            cublas_sync_stream();

            PROFSTART("Y:swap");

            k1 = jstart;
            k2 = MIN(m, jend );

            inc1 = descA[M_];
            inc2 = descA[M_];
            for(int kk=k1; kk <= k2; kk++) {

              ip = gipiv(kk);
              assert( (1 <= kk) && (kk <= ip) && (ip <= m) );

              if (ip != kk) {
                 /*
                  * interchange to left side A(kk,1:(jstart-1))
                  */
                 nn = (jstart-1);
                 i1 = (ia-1) + kk;
                 i2 = (ia-1) + ip;
                 j1 = (ja-1) + 1;
                 j2 = j1;
                 if (nn >= 1) {
                   scalapack_pcswap( &nn,  A,&i1,&j1,descA,&inc1,
                                           A,&i2,&j2,descA,&inc2);
                   };

                 /*
                  * interchange to right side A(ii,(jend+1):n
                  */
                 nn = n - (jend + 1) + 1;
                 i1 = (ia-1) + kk;
                 i2 = (ia-1) + ip;
                 j1 = (ja-1) + (jend+1);
                 j2 = j1;
                 if (nn >= 1) {
                   scalapack_pcswap( &nn, A,&i1,&j1,descA,&inc1,
                                          A,&i2,&j2,descA,&inc2 );
                 };
              };
            };
            PROFEND("Y:swap");

         }; /* end for jstart */

#ifdef USE_CUBLASV2
    {
      cudaError_t ierr;
      ierr = cudaFree( (void *) dY );
      assert( ierr == cudaSuccess );
      dY = 0;
    }
#else
         cu_status = cublasFree(dY);
         CHKERR(cu_status);
#endif


 /*
  * copy pivot vector 
  */
 if (ia != 1) {
  for(int i=1; i <= MIN(m,n); i++) {
   gipiv(i) = gipiv(i) + (ia-1);
   };
 };

 i1 = 1; 
 j1 = 1;
 inc1 = 1;

 i2 = (ia-1) + 1;
 j2 = (ja-1) + 1;
 inc2 = 1;

 lld = descA[MB_] + 
   Cnumroc( descA[M_], descA[MB_], myprow, descA[RSRC_], nprow );

 Cdescset( desc_ipiv, 
     descA[M_], 1, 
     descA[MB_], 1,
     descA[RSRC_], -1,
     descA[CTXT_], lld );

 PROFSTART("Y:ipiv");

 use_replicated_storage = FALSE;
 if (use_replicated_storage) {
   desc_gipiv[RSRC_] = -1;
   desc_gipiv[CSRC_] = -1;
   desc_ipiv[CSRC_] = -1;

   scalapack_picopy( &minmn, 
     &(gipiv(1)),&i1,&j1,desc_gipiv, &inc1,
     &(ipiv(1)),&i2,&j2,desc_ipiv, &inc2 );
   }
 else {
   int icontxt;
   int lrindx1,lcindx1,lrindx2,lcindx2,rsrc,csrc;
   char scope = 'R';
   char top = ' ';
   int msize,nsize,lld,lr1,lc1,lr2,lc2;

   icontxt = descA[CTXT_];
   rsrc = ia_proc;
   csrc = ja_proc;

   desc_gipiv[RSRC_] = rsrc;
   desc_gipiv[CSRC_] = csrc;
   desc_ipiv[CSRC_] = csrc;

   scalapack_picopy( &minmn, 
     &(gipiv(1)),&i1,&j1,desc_gipiv, &inc1,
     &(ipiv(1)),&i2,&j2,desc_ipiv, &inc2 );


   local_extent( minmn,n, ia,ja,descA,
          &msize,&nsize, &lr1,&lc1, &lr2,&lc2 );
   lrindx1 = lr1;
   lrindx2 = lr2;

   if (msize >= 1) {
    lld = msize;
    nsize = 1;
    if (npcol > 1) {
     if (mypcol == csrc) {
      scalapack_igebs2d(&icontxt, &scope, &top, 
         &msize, &nsize, &(ipiv(lrindx1)), &lld );
      }
     else {
      rsrc = myprow;
      scalapack_igebr2d(&icontxt, &scope, &top, 
         &msize, &nsize, &(ipiv(lrindx1)), &lld,  &rsrc,&csrc );
      };
    };
   };

 };
 PROFEND("Y:ipiv");


 if (gipiv_ != 0) {
   free(gipiv_); gipiv_ = 0;
 };

 if (ipivY_ != 0) {
   free(ipivY_); ipivY_ = 0;
 };

 if (use_new_hY) {
  if (hY != 0) {
   free(hY); hY = 0;
  };
 };
        
        return;
}

#ifdef __cplusplus
extern "C"
#endif
void pcgetrf_ooc2_( int *m_in, int *n_in,
                 float *A, int *ia_in, int *ja_in, int *descA, 
                          int *ipiv_,  int *memsize_in, int *info )
{
  pcgetrf_ooc2( m_in, n_in,
               A, ia_in, ja_in, descA,
               ipiv_, memsize_in, info );
}
