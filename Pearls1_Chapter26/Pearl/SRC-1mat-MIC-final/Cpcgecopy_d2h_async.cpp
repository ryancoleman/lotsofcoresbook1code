#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cpcgecopy_general_async(int m, int n, 
        void *A, int ia, int ja, int *descA,
        void *B, int ib, int jb, int *descB, int is_device_to_host)
{
#define dA(i,j)  (((cuComplex*)A) + IDX2F(i,j,descA[LLD_]))
#define dT(i,j) (((cuComplex *)T) + IDX2F(i,j,descT[LLD_]))

#define dB(i,j)  (((cuComplex *)B) + IDX2F(i,j,descB[LLD_]))
/*
  perform    copy

  B( ib:(ib+m-1), jb:(jb+n-1)) <-  A( ia:(ia+m-1),ja:(ja+n-1))

 */

const int use_MallocHost = FALSE;
const int use_igsum2d = FALSE;

cublasStatus cu_status;

cudaError_t cuda_status;

char notrans[] = "NoTrans";

int descT[DLEN_];

int ldA,ldB,ldT;

int is_same_context, is_same_mb, is_same_nb;
int is_same_p, is_same_q;
int is_same_offset;
int is_same_Locp, is_same_Locq;
int is_aligned;

int lrA1,lcA1, lrA2,lcA2;
int lrT1,lcT1, lrT2,lcT2;
int lrB1,lcB1, lrB2,lcB2;
int rsrc,csrc;
int rsrcA1,csrcA1,  rsrcA2, csrcA2;
int rsrcB1,csrcB1,  rsrcB2, csrcB2;
int iia,jja, iib,jjb;
int icontxt, nprow,npcol, myprow,mypcol;
int LocpA,LocqA,  LocpB,LocqB, LocpT,LocqT;
int mm,nn, lmm,lnn;
size_t nbytes;

float one_[REAL_PART+IMAG_PART+1];
float *one = &(one_[0]);

float zero_[REAL_PART+IMAG_PART+1];
float *zero = &(zero_[0]);


float alpha_[REAL_PART+IMAG_PART+1];
float *alpha = &(alpha_[0]);

float beta_[REAL_PART+IMAG_PART+1];
float *beta = &(beta_[0]);

int isize, isizeT;
float *T = 0;

int elemSize = sizeof(cuComplex);
int nnb, jstart,jend,jsize;
int is_ok;

int nmax;
const int bufsize =  1024*1024;
const int use_simple = FALSE;;

one[REAL_PART] = 1.0;
one[IMAG_PART] = 0.0;
zero[REAL_PART] = 0.0;
zero[IMAG_PART] = 0.0;

 if ((m <= 0) || (n <= 0)) {
   return;
 };

 T = 0; 

 ldA = descA[LLD_];
 ldB = descB[LLD_];

  icontxt = descA[CTXT_];
  Cblacs_gridinfo( icontxt, &nprow,&npcol, &myprow, &mypcol);
  assert( nprow >= 1);
  assert( npcol >= 1);
  assert( (0 <= myprow) && (myprow < nprow));
  assert( (0 <= mypcol) && (mypcol < npcol));

  is_ok = (1 <= ia) && (ia + m-1 <= descA[M_]);
  if (!is_ok) {
    printf("Cpcgecopy (%d,%d) :ia %d m %d descA[M_] %d  \n",
            myprow,mypcol,     ia,   m,   descA[M_] );
    printf("Cpcgecopy (%d,%d) :ja %d n %d descA[N_] %d \n",
            myprow,mypcol,     ja,   n,   descA[N_] );
    printf("Cpcgecopy (%d,%d) :ib %d jb %d descB[M_] %d descB[N_] %d\n",
            myprow,mypcol,     ib,   jb,   descB[M_],   descB[N_] );
  };
  assert( (1 <= ia) && (ia + m-1 <= descA[M_]));
  assert( (1 <= ja) && (ja + n-1 <= descA[N_]));
  assert( (1 <= ib) && (ib + m-1 <= descB[M_]));
  assert( (1 <= jb) && (jb + n-1 <= descB[N_]));


  is_same_context = (descA[CTXT_] == descB[CTXT_]);
  is_same_mb = (descA[MB_] == descB[MB_]);
  is_same_nb = (descA[NB_] == descB[NB_]);

  is_same_p = (Cindxg2p(ia,descA[MB_], myprow, descA[RSRC_],nprow) ==
               Cindxg2p(ib,descB[MB_], myprow, descB[RSRC_],nprow) );

  is_same_q = (Cindxg2p(ja,descA[NB_], mypcol, descA[CSRC_],npcol) ==
               Cindxg2p(jb,descB[NB_], mypcol, descB[CSRC_],npcol) );

  is_same_offset = (MOD(ia,descA[MB_]) == MOD(ib,descB[MB_])) &&
                   (MOD(ja,descA[NB_]) == MOD(jb,descB[NB_]));


  local_extent( m,n, ia,ja,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );


  local_extent( m,n, ib,jb,descB, &LocpB,&LocqB,&lrB1,&lcB1, &lrB2,&lcB2 );



  /*
  if ((LocpA >= 1) || (LocpB >= 1)) {
     is_same_Locp = (LocpA == LocpB);
  };
  if ((LocqA >= 1) || (LocqB >= 1)) {
     is_same_Locq = (LocqA == LocqB);
  };
  */

  is_same_Locq = (LocqA == LocqB);

  is_same_Locp = (LocpA == LocpB);

  is_aligned = is_same_context &&
               is_same_mb && is_same_nb &&
               is_same_p && is_same_q &&
               is_same_offset &&
               is_same_Locp && is_same_Locq;

  assert( is_same_q );

  assert( is_same_p );

  assert( is_same_offset );

  assert( is_same_Locp );

  assert( is_same_Locq );

  assert( is_aligned );


        
       /*
        no communication required
        copy from device to host
        */

       ldA = descA[LLD_];
       ldB = descB[LLD_];

       mm = LocpA;
       nn = LocqA;

       if (is_device_to_host) {
         /* 
          * transfer from device to host
          */
         if ( (mm >= 1) && (nn >= 1) ) {
#ifdef USE_CUBLASV2
           {
             cublasStatus_t istatus;
             istatus = cublasGetMatrixAsync(mm, nn, elemSize, 
                 (void *) dA(lrA1,lcA1), ldA, (void *) dB(lrB1,lcB1), ldB,
                 cublas_get_stream() );
             assert( istatus == CUBLAS_STATUS_SUCCESS );
           }
#else
           cu_status = cublasGetMatrix(mm,nn, elemSize,
               (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
            CHKERR(cu_status);
#endif
           };
         }
       else {
         /* 
          * transfer from host to device
          */
         if ( (mm >= 1) && (nn >= 1) ) {
#ifdef USE_CUBLASV2
           {
             cublasStatus_t istatus;

             istatus = cublasSetMatrixAsync(mm,nn,elemSize,
               (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB,
               cublas_get_stream()   );

             assert( istatus == CUBLAS_STATUS_SUCCESS );

           }
#else
            cu_status = cublasSetMatrix(mm,nn,elemSize,
               (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
            CHKERR(cu_status);
#endif
           };
         };
               



  return;
}

#ifdef __cplusplus
extern "C"
#endif
void Cpcgecopy_d2h_async(int m, int n, 
        cuComplex *A, int ia, int ja, int *descA,
        float *B, int ib, int jb, int *descB  )
{
  /*
   * Transfer device to host
   */

   int ltrue = (1 == 1);
   int is_device_to_host = ltrue;
   Cpcgecopy_general_async(m,n,A,ia,ja,descA,B,ib,jb,descB, is_device_to_host );
} 

#ifdef __cplusplus
extern "C"
#endif
void Cpcgecopy_h2d_async(int m, int n, 
        float *A, int ia, int ja, int *descA,
        cuComplex *B, int ib, int jb, int *descB  )
{
  /*
   * Transfer host to device
   */
   int lfalse = (1 == 0);
   int is_device_to_host = lfalse;
   Cpcgecopy_general_async(m,n,A,ia,ja,descA,B,ib,jb,descB, is_device_to_host );
} 
