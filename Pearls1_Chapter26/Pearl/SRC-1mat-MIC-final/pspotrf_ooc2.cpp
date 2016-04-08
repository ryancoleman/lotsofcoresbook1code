#include "ooclu.h"


#ifdef __cplusplus
extern "C" 
#endif
void pspotrf_ooc2( char *uplo_in, int *n_in, 
             float *A, int *ia_in, int *ja_in, int *descA, 
	     int *memsize_in, int *info )
{

int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;
int memsize = *memsize_in;

int use_replicated_storage = FALSE;
const int use_pspotrf_gpu2 = TRUE;

const int idebug = 1;
const int use_immediate_copy_to_host = TRUE;


char left[] = "Left";
char lower[] = "Lower";
char unit[] = "Unit";
char notrans[] = "NoTrans";

char *side = left;
char *uplo = lower; // pspotrf_ooc2 only work for lower triangular matrix now
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
float *dY = 0;

int mb,nb, nbx, nby, npanels, isize, isizeX,isizeY;
int iy0,jy0;
int ix0,jx0, iix,jjx, ix,iy, jx,jy;
int descY[DLEN_];

int descX_[DLEN_];
int *descX = &(descX_[0]);

int k1,k2,incx,ip;
int jstart,jend,jsize,i,j,jb, jstartm1;
int kfinal, kstart,kend,ksize, knb, ihY,jhY;
int isizehY;
int has_work;

int mm_lu, nn_lu,iy_lu,jy_lu;
int msize, minmn,mm, nn, kk, iia,jja, ii,jj, lld;
int iia_proc,jja_proc;
int iproc, rsrc, csrc, Locp, Locq, irsrc,icsrc;

int elemSize = sizeof(float);

float *hY = 0;
int deschY[DLEN_];
size_t nbytes;
int i1,j1,inc1,  i2,j2,inc2, lrindx;

const int use_new_hY = FALSE;
const int use_fixed_size_y_panel = FALSE;

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
 * perform Cholesky factorization similar to scalapack PSPOTRF
 * but use memsize entries on GPU
 */


/*
 * check arguments
 */

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


mm = n;
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
  printf("pspotrf_ooc2: mm %d nn %d mb %d nb %d isize %d \n",
                        mm,nn, mb,nb, isize );
  printf("pspotrf_ooc2: memsize %d npanels %d nby %d\n",
      memsize, npanels, nby );
};

minmn = n;

/*
 * directly use host matrix as "X" panel
 * no need to allocate storage for "X" panel
 */

 dX = A;
 ix = ia;
 jx = ja;
 descX = descA;


 if (use_fixed_size_y_panel) {
   /*
    * Preallocate storage for Y panel 
    */
   setup_desc( n,nby, ia,ja,descA,  &isizeY, descY );
   nbytes = isizeY;
   nbytes *= elemSize;
   }
  else {
    isizeY = memsize;
    nbytes = isizeY;
    nbytes *= elemSize;
    }


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


 /*
  * Main outer loop over panels
  */
 for(jstart=1; jstart <= n; jstart = jend + 1 ) {
    /*
     Adjust descY to be aligned with host data
     Y(1,1) is on same processor as A(iia,jja)
     */

      iia = (ia-1) + jstart;
      jja = (ja-1) + jstart;
      mm = n - jstart + 1;


    if (use_fixed_size_y_panel) {
        /* nby is unchanged, do nothing */
        }
    else {
 

      /*
       * For the same amount of device memory
       * try to use a wider Y-panel
       */
       
       /*
        * estimate the amount of memory needed for a mm by (nb * npcol) panel
        * and see how many copies of that will fit
        */
       nn = descA[NB_] * npcol;
       setup_desc(mm, nn, iia, jja, descA, &isize, descY );

       npanels = MAX(1, (int) (isizeY / isize ));
       nby = npanels * nn;
       };

       setup_desc( mm, nby, iia, jja, descA, &isize, descY );


    
    j = jstart;
    jend = MIN( n, jstart + nby - 1 );
    jsize = jend - jstart + 1;
    jb = jsize;

    /*
     Note index in Y panel starts at (1,1)
     
     Note we need to copy  only the lower part of matrix

     copy Y(1:isize,1:jsize) <-  A(jstart:n,jstart:jend)
     */

     iy0 = 1; 
     jy0 = 1;

     iy = (iy0-1) + 1;
     jy = (jy0-1) + 1;
     iia = (ia-1) + jstart;
     jja = (ja-1) + jstart;
     isize = n - jstart + 1;
     mm = isize;
     nn = jsize;

     PROFSTART("Y <- A");
     Cpsgecopy_h2d( mm,nn, A,iia,jja,descA,   dY, iy,jy,descY );
     PROFEND("Y <- A");


     /*
      Incorporate previous  updates into Y
      */
     kfinal = (ja-1) + (jstart-1);
     for(int kstart=ja; kstart <= kfinal;  kstart = kend + 1) {

       /*
        align kend on end of block boundary
        */
       kend = kstart + nbx - MOD(kstart-1,nbx) -1;
       kend = MIN(kend, kfinal);

       ksize = kend - kstart + 1;
       uplo = lower;
       trans = notrans;
       mm = isize;
       nn = jsize;
       kk = ksize;
       iia = (ia-1) + jstart;
       jja = (ja-1) + kstart;

       iy = (iy0-1) + 1;
       jy = (jy0-1) + 1;
       if ((mm >= 1) && (nn >= 1) && (kk >= 1)) {
         PROFSTART("pspotrf_ooc2:pdsyrk");

         Cpssyrk_hhd( *uplo, *trans, mm,nn,kk,
           neg_one,  A, iia,jja, descA, 
           one,      dY, iy, jy, descY );

         PROFEND("pspotrf_ooc2:pdsyrk");
         };
        };



/*
 *           -------------------------------------
 *           all left updates applied to panel d_Y
 *           perform factorization of  panel Y 
 *           in GPU device memory
 *           -------------------------------------
 */
            j = jstart;
            jsize = jend - jstart + 1;
            jb  = jsize;
            mm = n - j + 1;
            nn = jsize;

             /*
              * similar to pspotrf( &mm,&nn,Y,j,1,&descY, &iinfo);
             */
             iy = (iy0-1) + 1; 
             jy = (jy0-1) + 1;
             has_work = (mm >= 1) && (nn >= 1);
             if (has_work) {
               iinfo = 0;
               mm_lu = mm;
               nn_lu = nn;
               iy_lu = iy;
               jy_lu = jy;


                 int iah0 = (ia-1) + jstart;
                 int jah0 = (ja-1) + jstart;

                 PROFSTART("Y:pspotrf_gpu2");
                 pspotrf_gpu2( uplo, &mm_lu, &nn_lu, 
                   dY,&iy_lu,&jy_lu,descY, 
                   A, &iah0, &jah0, descA,
                   &iinfo );
                 PROFEND("Y:pspotrf_gpu2");
             };


#ifdef USE_FAKE_CUBLAS
             if (idebug >= 2) {
               char cmatnm[] = "dY2";
               if (is_root) { printf("after pspotrf_gpu\n"); };
               Cpslaprnt( mm,nn, (float *) dY,iy,jy,descY,cmatnm);
             };
#endif

             if (has_work) {
               if (iinfo < 0) {
                 *info = iinfo;
                 if (idebug >= 1) {
                   printf("pspotrf_ooc2: pspotrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
               if (iinfo > 0) {
                 *info = iinfo + (j-1);
                 if (idebug >= 1) {
                   printf("pspotrf_ooc2: pspotrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
             };


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


        return;
}

#ifdef __cplusplus
extern "C"
#endif
void pspotrf_ooc2_( char *uplo_in, int *n_in,
                 float *A, int *ia_in, int *ja_in, int *descA, 
                          int *memsize_in, int *info )
{
  pspotrf_ooc2( uplo_in, n_in,
               A, ia_in, ja_in, descA,
               memsize_in, info );
}
