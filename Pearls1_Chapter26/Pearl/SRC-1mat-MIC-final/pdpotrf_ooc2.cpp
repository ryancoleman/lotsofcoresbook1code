/*
Copyright (c) The University of Tennessee.  All rights reserved.


$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution.

- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against recipient by any third party for infringement of that parties intellectual property rights.
*/
#include "ooclu.h"


#ifdef __cplusplus
extern "C" 
#endif
void pdpotrf_ooc2( char *uplo_in, int *n_in, 
             double *A, int *ia_in, int *ja_in, int *descA, 
	     int *memsize_in, int *info )
{

int n = *n_in;
int ia = *ia_in;
int ja = *ja_in;
int memsize = *memsize_in;

int use_replicated_storage = FALSE;
const int use_pdpotrf_gpu2 = TRUE;
const int use_pdgemm_copy = FALSE;

const int idebug = 1;
const int use_immediate_copy_to_host = TRUE;


char left[] = "Left";
char lower[] = "Lower";
char unit[] = "Unit";
char transpose[] = "Transpose";
char ctranspose[] = "ConjuateTranspose";
char notrans[] = "NoTrans";

char *side = left;
char *uplo = lower; // pdpotrf_ooc2 only work for lower triangular matrix now
char *trans = notrans;
char *diag = unit;

char forward[] = "Forward";
char rowwise[] = "Rowwise";
char *direct = forward;
char *rowcol = rowwise;

#ifndef USE_MIC
cublasStatus cu_status;
#endif

int nprow,npcol,myprow,mypcol, icontxt;
int is_root;


double *X = 0; // triangular part of Y
double *dY = 0; // Y-panel on device

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

int elemSize = sizeof(double);

double *hY = 0;
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

double beta_[REAL_PART+IMAG_PART+1];
double *beta = &(beta_[0]);

double alpha_[REAL_PART+IMAG_PART+1];
double *alpha = &(alpha_[0]);

double one_[REAL_PART+IMAG_PART+1];
double *one = &(one_[0]);

double zero_[REAL_PART+IMAG_PART+1];
double *zero = &(zero_[0]);

double neg_one_[REAL_PART+IMAG_PART+1];
double *neg_one = &(neg_one_[0]);


zero[REAL_PART] = 0.0;
zero[IMAG_PART] = 0.0;

one[REAL_PART] = 1.0;
one[IMAG_PART] = 0.0;

neg_one[REAL_PART] = -1;
neg_one[IMAG_PART] = 0.0;


X = 0;
dY = 0;
hY = 0;

*info = 0;

/*
 * perform Cholesky factorization similar to scalapack PZPOTRF
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
  printf("pdpotrf_ooc2: mm %d nn %d mb %d nb %d isize %d \n",
                        mm,nn, mb,nb, isize );
  printf("pdpotrf_ooc2: memsize %d npanels %d nby %d\n",
      memsize, npanels, nby );
};

minmn = n;

/*
 * need to allocate storage for "X" panel
 */

 X = 0;

 if (use_fixed_size_y_panel) {
   /*
    * Preallocate storage for X panel 
    */
   setup_desc( nby,nby, ia,ja,descA,  &isizeX, descX );
   nbytes = isizeX;
   nbytes *= elemSize;
   }
 else {
   isizeX = memsize;
   nbytes = isizeX;
   nbytes *= elemSize;
 }

 X = (double*) malloc(nbytes);
 assert( X != 0 );

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

#ifdef USE_MIC
  dY = (double*) offload_Alloc(nbytes);
  assert( dY != 0 );
#else
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
        // isize is the same for all processes
       npanels = MAX(1, (int) (isizeY / isize ));
       nby = npanels * nn;
    };

    setup_desc( mm, nby, iia, jja, descA, &isize, descY );
    setup_desc( nby, nby, iia, jja, descA, &isizeX, descX );


   
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

     ix0 = 1;
     jx0 = 2;
     ix = iia;
     jx = jja + 1;
     iix = nn - 1;
     jjx = iix;

     PROFSTART("Y <- A");
       Cpdgecopy_h2d( mm,nn, A,iia,jja,descA,   dY, iy,jy,descY );
     PROFEND("Y <- A");
     PROFSTART("X <- A");
//   printf("X <- A skip no\n");
       scalapack_pdtradd( "U", "N", &iix, &jjx,
                          one, A, &ix, &jx, descA,
                          zero, X, &ix0, &jx0, descX);
     PROFEND("X <- A");


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
       ii = (ia-1) + kstart;
       jj = (ja-1) + jstart;

       iy = (iy0-1) + 1;
       jy = (jy0-1) + 1;
       if ((mm >= 1) && (nn >= 1) && (kk >= 1)) {

        if(use_pdgemm_copy){

         PROFSTART("pdpotrf_ooc2:pdgemm");
          PROFSTART("trans_copyA");
           trans = ctranspose;
           scalapack_pdgeadd( trans, &kk, &mm, 
              one, 
              A, &iia, &jja, descA,
              zero,
              A, &ii, &jj, descA );
          PROFEND("trans_copyA");

           trans = notrans;
           Cpdgemm_hhd( *trans, *trans, mm,nn,kk, 
               neg_one, A, iia,jja, descA, 
                  A, ii,jj,   descA, 
               one,  dY, iy,jy, descY );

         PROFEND("pdpotrf_ooc2:pdgemm");

        }else{

         PROFSTART("pdpotrf_ooc2:pdsyrk");

         Cpdsyrk_hhd( *uplo, *trans, mm,nn,kk,
           neg_one,  A, iia,jja, descA, 
           one,      dY, iy, jy, descY );

         PROFEND("pdpotrf_ooc2:pdsyrk");
         };
        };
     }



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
              * similar to pdpotrf( &mm,&nn,Y,j,1,&descY, &iinfo);
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

                 PROFSTART("Y:pdpotrf_gpu2");
                 pdpotrf_gpu2( uplo, &mm_lu, &nn_lu, 
                   dY,&iy_lu,&jy_lu,descY, 
                   A, &iah0, &jah0, descA,
                   &iinfo );
                 PROFEND("Y:pdpotrf_gpu2");
             };


#ifdef USE_FAKE_CUBLAS
             if (idebug >= 2) {
               char cmatnm[] = "dY2";
               if (is_root) { printf("after pdpotrf_gpu\n"); };
               Cpdlaprnt( mm,nn, (double *) dY,iy,jy,descY,cmatnm);
             };
#endif

             if (has_work) {
               if (iinfo < 0) {
                 *info = iinfo;
                 if (idebug >= 1) {
                   printf("pdpotrf_ooc2: pdpotrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
               if (iinfo > 0) {
                 *info = iinfo + (j-1);
                 if (idebug >= 1) {
                   printf("pdpotrf_ooc2: pdpotrf return iinfo=%d\n", iinfo);
                 };
                 return;
                 };
             };

     PROFSTART("A <- X");
       scalapack_pdtradd( "U", "N", &iix, &jjx,
                          one, X, &ix0, &jx0, descX,
                          zero, A, &ix, &jx, descA);
     PROFEND("A <- X");

     }; /* end for jstart */

#ifdef USE_MIC
    offload_Free(dY);
#else
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
#endif

        return;
}

#ifdef __cplusplus
extern "C"
#endif
void pdpotrf_ooc2_( char *uplo_in, int *n_in,
                 double *A, int *ia_in, int *ja_in, int *descA, 
                          int *memsize_in, int *info )
{
  pdpotrf_ooc2( uplo_in, n_in,
               A, ia_in, ja_in, descA,
               memsize_in, info );
}
