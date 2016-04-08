#include "scalapack.h"




#ifdef __cplusplus
extern "C" 
#endif
void pzgecopy_hd( F_CHAR_T TRANS, int *m_in, int *n_in,
    double *A, int *ia_in, int *ja_in, int *descA,
    double *dB, int *ib_in, int *jb_in, int *descB )
{
/*
 Copy m by n distributed submatrix  from host to GPU
 */

int m  = *m_in;
int n  = *n_in;

int ia = *ia_in;
int ja = *ja_in;
int ib = *ib_in;
int jb = *jb_in;

int nprow = 0;
int npcol = 0;
int myprow = 0;
int mypcol = 0;
int mmb = 0;
int nnb = 0;

int istart = 0;
int iend = 0;
int isize = 0;

int Locp = 0;
int Locq = 0;
int lld = 0;

int jstart = 0;
int jend = 0;
int jsize = 0;

int iib = 0;
int jjb = 0;

cuDoubleComplex *dBptr = 0;
double *Btmp = 0;
F_CHAR_T NoTrans = "N";

int descBtmp[DLEN_];
int elmSize = sizeof(cuDoubleComplex);

double z_one[2];
double z_zero[2];

/*
 Tuneable parameters
 */
int mfactor = 4;
int nfactor = 4;


int rsrc = 0;
int csrc = 0;
int irsrc = 0;
int jcsrc = 0;
int info = 0;

int notran = 0;
int TrA = 0;

int lrindx = 0;
int lcindx = 0;
int nrow = 0;
int ncol = 0;
int mm = 0;
int nn = 0;
int iia = 0;
int jja = 0;

cublasStatus cu_status;

z_one[REAL_PART] = 1;
z_one[IMAG_PART] = 0;
z_zero[REAL_PART] = 0;
z_zero[IMAG_PART] = 0;

Cblacs_gridinfo( descA[CTXT_], &nprow, &npcol, &myprow, &mypcol );

notran = ( ( TrA = Mupcase( F2C_CHAR( TRANS )[0] ) ) == CNOTRAN );


/*
 * check arguments
 */
if ((m <= 0) || (n <= 0)) {
   return;
   };

assert( (1 <= ia) && ((ia + m-1) <= descA[M_] ) );
assert( (1 <= ja) && ((ja + n-1) <= descA[N_] ) );
assert( (1 <= ib) && ((ib + m-1) <= descB[M_] ) );
assert( (1 <= jb) && ((jb + n-1) <= descB[N_] ) );

/*
 * Create a temp matrix that is aligned to descB.
 * Assume size is   mmb by nnb
 */

if (notran) {
  mmb = MIN( m, descB[MB_] * nprow * mfactor);
  nnb = MIN( n, descB[NB_] * npcol * nfactor);
  }
else {
  nnb = MIN( m, descB[MB_] * nprow * mfactor);
  mmb = MIN( n, descB[NB_] * npcol * nfactor);
  };
  

mmb = MAX( mmb,1);
nnb = MAX( nnb,1);

rsrc = indxg2p_( &ib, &descB[MB_], &myprow, &descB[RSRC_], &nprow);
csrc = indxg2p_( &jb, &descB[NB_], &mypcol, &descB[CSRC_], &npcol);


Locp = numroc_( &mmb, &descB[MB_], &myprow, &rsrc, &nprow );
Locq = numroc_( &nnb, &descB[NB_], &mypcol, &csrc, &npcol );

Btmp = (double*) malloc( MAX(1,(Locp * Locq))*elmSize );
assert( Btmp != 0 );


lld = MAX(1, Locp);
descinit_( descBtmp, &mmb, &nnb, &descB[MB_], &descB[NB_],
            &rsrc, &csrc, &descB[CTXT_], &lld, &info );

assert( info == 0);



for( jstart=ja; jstart <= ja + n-1; jstart = jend + 1) {
  jend = MIN( ja + n -1, jstart + nnb - 1);
  jsize = jend - jstart + 1;

  for( istart=ia; istart <= ia + m-1; istart = iend + 1) {
     iend = MIN( ia + m-1, istart + mmb -1);
     isize = iend - istart + 1;

     iia = ia + (istart-1);
     jja = ja + (jstart-1);


     iib = 1;
     jjb = 1;
     

     if (notran) {
        mm = isize;
        nn = jsize;
        }
    else {
        mm = jsize;
        nn = isize;
        };


     pzgeadd_( TRAN,  &mm, &nn,   z_one, A, &iia, &jja, descA,  
                         z_zero, Btmp, &iib, &jjb, descBtmp );

     /* 
      * find local extent
      */

     
     if (notran) {
       iib = ib + (istart-1);
       jjb = jb + (jstart-1);
       }
     else {
       iib = ib + (jstart-1);
       jjb = jb + (istart-1);
       };

     if (notran) {
       nrow = numroc_( &isize, &descB[MB_], &myprow, &rsrc, &nprow);
       ncol = numroc_( &jsize, &descB[NB_], &mypcol, &csrc, &npcol);
       }
     else {
       nrow = numroc_( &jsize, &descB[MB_], &myprow, &rsrc, &nprow);
       ncol = numroc_( &isize, &descB[NB_], &mypcol, &csrc, &npcol);

     };

     /*
      Perform global
      dB( iib:(iib+isize-1), jjb:(jjb+jsize-1)) <- B(1:isize,1:jsize)

      Perform local
      dB( lrindx:(lrindx+nrow-1), lcindx:(lcindx+ncol-1)) <-
              B(1:nrow, 1:ncol)
      */


     infog2l_( &iib, &jjb, descB, &nprow, &npcol, &myprow, &mypcol,
               &lrindx, &lcindx,  &irsrc, &jcsrc );


     dBptr = (cuDoubleComplex *) dB;
     dBptr = dBptr + INDX2F( lrindx,lcindx, descB[LLD_]);


     cu_status = cublasSetMatrix( nrow,ncol,elmSize,
               (cuDoubleComplex *) Btmp,  descBtmp[LLD_],
                                   dBptr, descB[LLD_] );
     assert( cu_status == CUBLAS_STATUS_SUCCESS );

     };
  };

  free( Btmp );

  return;
}

     
  

