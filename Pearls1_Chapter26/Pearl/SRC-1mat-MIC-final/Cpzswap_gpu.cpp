#include "ooclu.h"

#define dA(i,j)   (  ((cuDoubleComplex *)A) + IDX2F((i),(j),descA[LLD_]))
#define dB(i,j)   (  ((cuDoubleComplex *)B) + IDX2F((i),(j),descB[LLD_]))

#ifdef __cplusplus
extern "C"
#endif
void Cpzswap_gpu( int n, cuDoubleComplex *A, int ia,int ja,int *descA, int incA,
                    cuDoubleComplex *B, int ib,int jb,int *descB, int incB )
{
/*
 perform pzswap operation when
 both distributed arrays A and B are in device memory
 */



/*
 * allocate temporary space on host
 * then use pzswap for communication
 */

const int use_MallocHost = FALSE;

cublasStatus cu_status;
size_t nbytes;
int elemSize = sizeof( cuDoubleComplex );

double *Atmp = 0;
double *Btmp = 0;

int descAtmp[DLEN_];
int descBtmp[DLEN_];
int ldA, ldB, ldAtmp, ldBtmp;

int nprow,npcol,myprow,mypcol;
int Locp, Locq, lrindx, lcindx, mm,nn;
int LocpA, LocqA, lrindxA, lcindxA;
int LocpB, LocqB, lrindxB, lcindxB;
int isizeA, isizeB, rsrc, csrc;

int iia,jja, iib, jjb;
int incAtmp, incBtmp;
int lrA1,lcA1, lrA2,lcA2;
int lrB1,lcB1, lrB2,lcB2;

Cblacs_gridinfo( descA[CTXT_], &nprow, &npcol, &myprow, &mypcol );


/*
 * allocate storage for vector from A
 */

if (incA == 1) {
   /*
    *  This is a column vector
    */
   mm = n; nn = 1;
   }
else {
  /*
   * This is a row vector
   */
   mm = 1; nn = n;
   };
setup_desc(  mm,nn, ia,ja, descA,   &isizeA, descAtmp );

nbytes = elemSize;
nbytes *= isizeA;
if (use_MallocHost) {
  Atmp = (double *) MallocHost( nbytes );
  }
else {
  Atmp = (double *) malloc( nbytes );
  };
assert( Atmp != 0 );


/*
 * copy vector from A
 */

PROFSTART("swap:GetMatrix");

local_extent( mm,nn,ia,ja,descA,  &LocpA, &LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );

lrindxA = lrA1;
lcindxA = lcA1;

ldA = descA[LLD_];
ldAtmp = descAtmp[LLD_];
if ( (LocpA >= 1) && (LocqA >= 1)) {
  /*
   * copy from GPU device to host CPU
   */
  cu_status = cublasGetMatrix( LocpA,LocqA, elemSize,
             dA(lrindxA,lcindxA), ldA,  Atmp, ldAtmp );

  CHKERR(cu_status);
  };



/*
 * allocate storage for vector from B
 */

Cblacs_gridinfo( descB[CTXT_], &nprow, &npcol, &myprow, &mypcol );

if (incB == 1) {
   /*
    *  This is a column vector
    */
   mm = n; nn = 1;
   }
else {
  /*
   * This is a row vector
   */
   mm = 1; nn = n;
   };
setup_desc( mm,nn, ib,jb,descB, &isizeB, descBtmp );

ldBtmp = descBtmp[LLD_];
ldB = descB[LLD_];

nbytes = elemSize;
nbytes *= isizeB;
if (use_MallocHost) {
  Btmp = (double *) MallocHost( nbytes );
  }
else {
  Btmp = (double *) malloc( nbytes );
  };
assert( Btmp != 0 );



/*
 * copy vector from B
 */

local_extent( mm,nn,ib,jb,descB,  &LocpB, &LocqB, &lrB1,&lcB1,  &lrB2,&lcB2 );

lrindxB = lrB1;
lcindxB = lcB1;



ldB = descB[LLD_];
ldBtmp = descBtmp[LLD_];
if ((LocpB >= 1) && (LocqB >= 1)) {
  /*
   * Copy from GPU to CPU host
   */
  cu_status = cublasGetMatrix(LocpB,LocqB,elemSize,
         dB(lrindxB,lcindxB), ldB, Btmp, ldBtmp );
  CHKERR(cu_status );
  };

PROFEND("swap:GetMatrix");

iia = 1; jja = 1;
iib = 1; jjb = 1;
if (incA == 1) {
   incAtmp = 1;
   }
else {
  incAtmp = descAtmp[M_];
};

if (incB == 1) {
   incBtmp = 1;
    }
else {
   incBtmp = descBtmp[M_];
};


PROFSTART("swap:pzswap");
scalapack_pzswap( &n, Atmp, &iia, &jja, descAtmp, &incAtmp,
                      Btmp, &iib, &jjb, descBtmp, &incBtmp );
PROFEND("swap:pzswap");


/*
 * copy from host CPU back to GPU
 */

PROFSTART("swap:SetMatrix");

if ((LocpA >= 1) && (LocqA >= 1)) {
  /*
   * Copy from CPU host to GPU device
   */
  cu_status = cublasSetMatrix( LocpA, LocqA, elemSize,
              Atmp, ldAtmp, dA(lrindxA,lcindxA), ldA );
  CHKERR(cu_status);
  };


if ((LocpB >= 1) && (LocqB >= 1)) {
  /*
   * Copy from CPU host to GPU device
   */
  cu_status = cublasSetMatrix( LocpB, LocqB, elemSize,
                 Btmp, ldBtmp, dB(lrindxB,lcindxB), ldB );
  CHKERR(cu_status);
  };

PROFEND("swap:SetMatrix");

/*
 * clean up
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

if (Btmp != 0) {
  if (use_MallocHost) {
    FreeHost(Btmp);
     }
  else {
    free(Btmp); 
  };
  Btmp = 0;
 };



return;


}

