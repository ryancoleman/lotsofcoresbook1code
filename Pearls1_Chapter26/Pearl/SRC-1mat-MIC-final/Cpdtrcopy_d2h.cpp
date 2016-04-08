#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cpdtrcopy_general(char uplo, int m, int n, 
        void *A, int ia, int ja, int *descA,
        void *B, int ib, int jb, int *descB, int is_device_to_host)
{
#define dA(i,j)  (((double*)A) + IDX2F(i,j,descA[LLD_]))
#define dT(i,j) (((double *)T) + IDX2F(i,j,descT[LLD_]))

#define dB(i,j)  (((double *)B) + IDX2F(i,j,descB[LLD_]))
/*
  perform    copy

  B( ib:(ib+m-1), jb:(jb+n-1)) <-  A( ia:(ia+m-1),ja:(ja+n-1))
  the global matrices A,B are trapezoidal
 */

const int use_MallocHost = FALSE;
const int use_igsum2d = FALSE;

#ifndef USE_MIC
cublasStatus cu_status;
cudaError_t cuda_status;
#endif

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
int iia,jja, iib,jjb, iit, jjt;
int icontxt, nprow,npcol, myprow,mypcol;
int LocpA,LocqA,  LocpB,LocqB, LocpT,LocqT;
int mm,nn, lmm,lnn;
size_t nbytes;

double one_[REAL_PART+IMAG_PART+1];
double *one = &(one_[0]);

double zero_[REAL_PART+IMAG_PART+1];
double *zero = &(zero_[0]);


double alpha_[REAL_PART+IMAG_PART+1];
double *alpha = &(alpha_[0]);

double beta_[REAL_PART+IMAG_PART+1];
double *beta = &(beta_[0]);

int is_lower = ( uplo == 'L') || (uplo == 'l');
int is_upper = ( uplo == 'U') || (uplo == 'u');
    assert( is_lower ||  is_upper );
int has_diag;

int isize, isizeT;
double *T = 0;

int elemSize = sizeof(double);
int nnb, jstart,jend,jsize,jfinal;
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

  is_ok = (1 <= ia) && (ia + m-1 <= descA[M_]) &&
          (1 <= ja) && (ja + n-1 <= descA[N_]) &&
          (1 <= ib) && (ib + m-1 <= descB[M_]) &&
          (1 <= jb) && (jb + n-1 <= descB[N_]);

  if (!is_ok) {
    printf("Cpdtrcopy (%d,%d) :ia %d m %d descA[M_] %d  \n",
            myprow,mypcol,     ia,   m,   descA[M_] );
    printf("Cpdtrcopy (%d,%d) :ja %d n %d descA[N_] %d \n",
            myprow,mypcol,     ja,   n,   descA[N_] );
    printf("Cpdtrcopy (%d,%d) :ib %d jb %d descB[M_] %d descB[N_] %d\n",
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



  if ((LocpA >= 1) || (LocpB >= 1)) {
     is_same_Locp = (LocpA == LocpB);
  };
  if ((LocqA >= 1) || (LocqB >= 1)) {
     is_same_Locq = (LocqB == LocqB);
  };

  is_aligned = is_same_context &&
               is_same_mb && is_same_nb &&
               is_same_p && is_same_q &&
               is_same_offset &&
               is_same_Locp && is_same_Locq;
  /*
   * check that all processors agree that it is aligned
   */

  if (use_igsum2d) {
  char scope = 'A';
  char top = ' ';
  int naligned = 0;
  int mm = 1;
  int nn = 1;
  int lld = 1;
  int rdest = -1;
  int cdest = -1;

  naligned = 0;
  if (is_aligned) { naligned = 1; };

  scalapack_igsum2d( &icontxt, &scope, &top,
      &mm, &nn, &naligned, &lld, &rdest, &cdest );

  is_aligned = is_aligned && (naligned == (nprow*npcol));
  } else {
    /* 
     * just to be safe
       is_aligned = FALSE;
     */
  };



  if (is_aligned) {
        
       /*
        no communication required
        copy from device to host
        */

    ldA = descA[LLD_];
    ldB = descB[LLD_];

    if(is_lower){ 
        if(ja < ia){
            local_extent( m,MIN(ia-ja,n), ia,ja,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent( m,MIN(ia-ja,n), ib,jb,descB, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            n = n - MIN(ia-ja,n);
            mm = LocpA;
            nn = LocqA;
        }else{
            m = MAX(m - ja + ia,0);
            ia = ja; ib = jb;
            mm = 0;
            nn = 0;
        }
        n = MIN(n,m);
    }else{
        if(ia < ja){ 
            local_extent( MIN(ja-ia,m),n, ia,ja,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent( MIN(ja-ia,m),n, ib,jb,descB, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            m = m - MIN(ja-ia,m);
            mm = LocpA;
            nn = LocqA;
        }else{
            n = MAX(n - ia + ja,0);
            ja = ia; jb = ib;
            mm = 0;
            nn = 0;
        }
        m = MIN(m,n);
    }

    if (is_device_to_host) {
        /* 
         * transfer from device to host
         */
        if ( (mm >= 1) && (nn >= 1) ) {
        #ifdef USE_MIC
            offload_dGetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
        #else
            cu_status = cublasGetMatrix(mm,nn, elemSize,
                (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                CHKERR(cu_status);
        #endif
        };
    }else{
        /* 
         * transfer from host to device
         */
        if ( (mm >= 1) && (nn >= 1) ) {
        #ifdef USE_MIC
            offload_dSetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
        #else
            cu_status = cublasSetMatrix(mm,nn, elemSize,
                (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                CHKERR(cu_status);
        #endif
        };
    };

    for(jstart=1; jstart <= n; jstart = jend + 1) {
        iia = ia + jstart - 1; 
        jend = MIN(((iia-1)/descA[NB_] + 1)*descA[NB_],((iia-1)/descA[MB_] + 1)*descA[MB_]);
        jend = jend - iia + 1;
        jend = MIN(jend,n);
        jsize = jend-jstart+1;
 
        has_diag = ( Cindxg2p(iia, descA[MB_], myprow, descA[RSRC_], nprow) == myprow ) &&
                   ( Cindxg2p(iia, descA[NB_], mypcol, descA[CSRC_], npcol) == mypcol );

        if(is_lower){
            local_extent(m-jend,jsize, iia+jsize,iia,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent(m-jend,jsize, iia+jsize,iia,descB, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            if(has_diag){
                mm = LocpA + jsize;
                nn = LocqA;
                lrA1 = lrA1 - jsize;
                lrB1 = lrB1 - jsize;
            }else{
                mm = LocpA;
                nn = LocqA;
            }
        }else{
            local_extent(jsize,n-jend, iia,iia+jsize,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent(jsize,n-jend, iia,iia+jsize,descA, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            if(has_diag){
                mm = LocpA;
                nn = LocqA + jsize;
                lcA1 = lcA1 - jsize;
                lcB1 = lcB1 - jsize;
            }else{
                mm = LocpA;
                nn = LocqA;
            }
        }

        if ( (mm >= 1) && (nn >= 1) ) {
            if(has_diag){
                if (is_device_to_host) {
                    /* 
                     * transfer from device to host
                     */
                #ifdef USE_MIC
                    offload_dtrGetMatrix(uplo, mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
                #else
                    cu_status = cublasGetMatrix(mm,nn, elemSize,
                        (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                        CHKERR(cu_status);
                #endif
                }else{
                    /* 
                     * transfer from host to device
                     */
                #ifdef USE_MIC
                    offload_dtrSetMatrix(uplo, mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
                #else
                    cu_status = cublasSetMatrix(mm,nn, elemSize,
                        (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                        CHKERR(cu_status);
                #endif
                }
            }else{
                if (is_device_to_host) {
                    /* 
                     * transfer from device to host
                     */
                #ifdef USE_MIC
                    offload_dGetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
                #else
                    cu_status = cublasGetMatrix(mm,nn, elemSize,
                        (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                        CHKERR(cu_status);
                #endif
                }else{
                    /* 
                     * transfer from host to device
                     */
                #ifdef USE_MIC
                    offload_dSetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
                #else
                    cu_status = cublasSetMatrix(mm,nn, elemSize,
                        (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                        CHKERR(cu_status);
                #endif
                }
            }
        }
    }
  } else {
    /*
     * Need communication or more complicated handling
     */

    /*
     * Try to limit transfer to approximately "bufsize" 
     * entries on each local processor
     */
    
    if(is_device_to_host){


    /* 
     * wernghurenhuslngvhusldnhurgsrgounvshlrgnvh
     * aeiognreiuauiaeonrghupanehrugpnaehrgohdjfkglwergbfsmdgshjgwslerguslerhguslrg
     */



    }else{
      if(is_lower){ 
        nnb = descB[NB_]*npcol;
        setup_desc(mm, nnb, 1, 1, descB, &isize, descT);
      }else{
        nnb = descB[MB_]*nprow;
        setup_desc(nnb, n, 1, 1, descB, &isize, descT);
      }
      isizeT = isize;
      nbytes = isizeT;
      nbytes *= elemSize;
      dT = (double*) malloc(nbytes);

      ldA = descA[LLD_];
      ldB = descB[LLD_];
      ldT = descT[LLD_];

      if(is_lower){
        if(ja < ia){
            jfinal = MIN(ia-ja,n);
        }else{
            jfinal = 0;
            ib = ib - (ja - ia);
            ia = ja;
        } 
      }else{
        if(ia < ja){
            jfinal = MIN(ja-ia,m);
        }else{
            jfinal = 0;
            jb = jb - (ia - ja);
            ja = ia;
        }
      }

      for(jstart=1; jstart <= jfinal; jstart = jend + 1) {
          if(is_lower){
            iib = ib;
            jjb = jb + jstart - 1;
            jend = MIN(((jjb-1)/nnb + 1)*nnb,nn);
            jsize = jend - jstart + 1; 
            iia = ia;
            jja = ja + jstart - 1;
            iit = 1;
            jjt = (jjb-1)%nnb + 1;

            mm = m;
            nn = jsize;
          }else{
            iib = ib + jstart - 1;
            jjb = jb;
            jend = MIN(((iib-1)/nnb + 1)*nnb,nn);
            jsize = jend - jstart + 1; 
            iia = ia + jstart - 1;
            jja = ja;
            iit = (iib-1)%nnb + 1;
            jjt = 1;

            mm = jsize;
            nn = n;
          }


          scalapack_pdgeadd(notrans, &mm, &nn,
              one, dA, &iia, &jja, descA,
              zero, dT, &iit, &jjt, descT);

          local_extent(mm,nn, iit,jjt,descT, &LocpT,&LocqT, &lrT1,&lcT1, &lrT2,&lcT2 );
          local_extent(mm,nn, iib,jjb,descB, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
          
          assert( LocpT == LocpB );
          assert( LocqT == LocqB );

          lmm = LocpT;
          lnn = LocqT;

          if ((lmm >= 1) && (lnn >= 1)) {
          #ifdef USE_MIC
            offload_dSetMatrix(lmm, lnn,
                dT(lrT1,lcT1),ldT, dB(lrB1,lcB1), ldB);
          #else
            cu_status = cublasSetMatrix(lmm,lnn,elemSize,
                dT(lrT1,lcT1),ldT, dB(lrB1,lcB1), ldB);

            is_ok = (cu_status == CUBLAS_STATUS_SUCCESS);
            if (!is_ok) {
                printf("Cpdtrcopy(%d,%d): m %d n %d lmm %d lnn %d isizeT %d lrT1 %d lcT1 %d lrB1 %d lcB1 %d iia %d jja %d  iib %d jjb %d \n",
                    myprow, mypcol, m, n,  lmm,lnn,
                    isizeT, lrT1,lcT1,   lrB1,lcB1,   iia,jja,  iib,jjb );
            };
            CHKERR(cu_status);
          #endif
          };
      }

      if(is_lower){
        n = n - jfinal;
        ja = ja + jfinal;
        jb = jb + jfinal;
      }else{
        m = m - jfinal;
        ia = ia + jfinal;
        ib = ib + jfinal;
      }

      for(jstart=1; jstart <= n; jstart = jend + 1) {
        iib = ib;
        jjb = jb + jstart - 1;
        jend = MIN(((iia-1)/descA[NB_] + 1)*descA[NB_],((iia-1)/descA[MB_] + 1)*descA[MB_]);
        jend = jend - iia + 1;
        jend = MIN(jend,n);
        jsize = jend-jstart+1;
        iia = ia;
        jja = ja + jstart - 1;
        iit = 1;
        jjt = (jjb-1)%nnb + 1;
 
        has_diag = ( Cindxg2p(iia, descA[MB_], myprow, descA[RSRC_], nprow) == myprow ) &&
                   ( Cindxg2p(iia, descA[NB_], mypcol, descA[CSRC_], npcol) == mypcol );

        if(is_lower){
            local_extent(m-jend,jsize, iia+jsize,iia,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent(m-jend,jsize, iia+jsize,iia,descB, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            if(has_diag){
                mm = LocpA + jsize;
                nn = LocqA;
                lrA1 = lrA1 - jsize;
                lrB1 = lrB1 - jsize;
            }else{
                mm = LocpA;
                nn = LocqA;
            }
        }else{
            local_extent(jsize,n-jend, iia,iia+jsize,descA, &LocpA,&LocqA, &lrA1,&lcA1, &lrA2,&lcA2 );
            local_extent(jsize,n-jend, iia,iia+jsize,descA, &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2 );
            if(has_diag){
                mm = LocpA;
                nn = LocqA + jsize;
                lcA1 = lcA1 - jsize;
                lcB1 = lcB1 - jsize;
            }else{
                mm = LocpA;
                nn = LocqA;
            }
        }

        if ( (mm >= 1) && (nn >= 1) ) {
            if(has_diag){
                /* 
                 * transfer from host to device
                 */
            #ifdef USE_MIC
                offload_dSetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
            #else
                cu_status = cublasSetMatrix(mm,nn, elemSize,
                    (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                    CHKERR(cu_status);
            #endif
            }else{
                /* 
                 * transfer from host to device
                 */
            #ifdef USE_MIC
                offload_dSetMatrix(mm, nn, dA(lrA1,lcA1), ldA, dB(lrB1,lcB1), ldB);
            #else
                cu_status = cublasSetMatrix(mm,nn, elemSize,
                    (void *) dA(lrA1,lcA1), ldA,  (void *) dB(lrB1,lcB1),ldB );
                    CHKERR(cu_status);
            #endif
            }
        }
      }






    }

    if (is_device_to_host) {
      int isize, npanel, desc_dummy[DLEN_];

      if (use_simple) {
         nnb = MIN(n,npcol*descA[NB_]);
         if (m == 1) { nnb = n; };

      } else {
         setup_desc( m, npcol*descA[NB_], ia,ja,descA, &isize, desc_dummy );
         npanel = MAX(1, bufsize/isize );
         nmax = npanel*(npcol*descA[NB_]);
         nnb = MIN(n,nmax);
      };
    } else {
      /*
       * Host to device
       */
      int isize, npanel,desc_dummy[DLEN_];

      if (use_simple) {
        nnb = MIN(n, npcol*descB[NB_]);
        if (m == 1) { nnb = n; };
      } else {
         setup_desc(m, npcol*descB[NB_], ib,jb,descB, &isize, desc_dummy);
         npanel = MAX(1, bufsize/isize);
         nmax = npanel*(npcol*descB[NB_]);
         nnb = MIN(n,nmax);
      };
    };


    /* 
     * preallocate T
     */

    if (is_device_to_host) {
       setup_desc(m,nnb, ia,ja,descA, &isize,descT );
    } else {
      setup_desc(m,nnb, ib,jb,descB, &isize, descT );
    };
    isizeT = isize;

    nbytes = isizeT;
    nbytes *= elemSize;
    T = (double *) MallocHost( nbytes );
    assert( T != 0 );





    for(jstart=1; jstart <= n; jstart = jend + 1) {
         jend = MIN(n, jstart + nnb-1);
         jsize = jend - jstart + 1;


         mm = m;
         nn = jsize;

         if (is_device_to_host) {
           iia = (ia-1) + 1;
           jja = (ja-1) + jstart;
           setup_desc( mm,nn, iia,jja, descA,  &isize, descT );
           }
         else {
           iib = (ib-1) + 1;
           jjb = (jb-1) + jstart;
           setup_desc( mm,nn, iib,jjb, descB, &isize, descT );
          };

         assert(isizeT >=  isize);


         ldA = descA[LLD_];
         ldB = descB[LLD_];
         ldT = descT[LLD_];

         iia = (ia-1) + 1;
         jja = (ja-1) + jstart;
         local_extent(mm,nn,iia,jja,descA, 
             &LocpA,&LocqA,&lrA1,&lcA1, &lrA2,&lcA2);

         iib = (ib-1) + 1;
         jjb = (jb-1) + jstart;
         local_extent(mm,nn,iib,jjb,descB, 
             &LocpB,&LocqB, &lrB1,&lcB1, &lrB2,&lcB2);

         local_extent(mm,nn,1,1,descT,  
                      &LocpT, &LocqT,
                      &lrT1,&lcT1, &lrT2,&lcT2);


      
         if (is_device_to_host) {

           assert( LocpT == LocpA );
           assert( LocqT == LocqA );

           lmm = LocpT;
           lnn = LocqT;;

           /*
            * transfer from A on GPU to algined matrix T on CPU
            */
           if ((lmm >= 1) && (lnn >= 1)) {
           #ifdef USE_MIC
            offload_dGetMatrix(lmm, lnn,
                dA(lrA1,lcA1), ldA, dT(lrT1,lcT1), ldT);
           #else
            cu_status = cublasGetMatrix(lmm,lnn,elemSize,
               dA(lrA1,lcA1), ldA,  dT(lrT1,lcT1),ldT );
            CHKERR(cu_status);
           #endif
             };
        
         iia = 1; 
         jja = 1;
         iib = (ib-1) + 1; 
         jjb = (jb-1) + jstart;

         
         alpha = one;
         beta = zero;


         if ((mm >= 1) && (nn >= 1)) {
           scalapack_pdgeadd( notrans, &mm, &nn, 
              alpha,  
              T, &iia,&jja, descT,
              beta,   
              (double *) B,    &iib,&jjb, descB );
            };
          }
         else {

            /*
             * transfer from host to device
             */


            /*
             * copy to matrix T, aligned to GPU data
             */
            alpha = one;
            beta = zero;

            iia = ia;
            jja = (ja-1) + jstart;
            iib = 1;
            jjb = 1;
            if ((mm >= 1) && (nn >= 1)) {
              scalapack_pdgeadd( notrans, &mm, &nn,
                alpha,       (double *) A,&iia,&jja, descA,
                beta,        T,&iib,&jjb, descT );
            };
                 

            assert( LocpT == LocpB );
            assert( LocqT == LocqB );

            /*
             * transfer from T on host to B on GPU
             */
            lmm = LocpT; lnn = LocqT;
            if ((lmm >= 1) && (lnn >= 1)) {
            #ifdef USE_MIC
               offload_dSetMatrix(lmm, lnn,
                    dT(lrT1,lcT1),ldT, dB(lrB1,lcB1), ldB);
            #else
               cu_status = cublasSetMatrix(lmm,lnn,elemSize,
                    dT(lrT1,lcT1),ldT,   dB(lrB1,lcB1), ldB );

               is_ok = (cu_status == CUBLAS_STATUS_SUCCESS);
               if (!is_ok) {
                 printf("Cpdgecopy(%d,%d): m %d n %d lmm %d lnn %d isizeT %d lrT1 %d lcT1 %d lrB1 %d lcB1 %d iia %d jja %d  iib %d jjb %d \n",
                     myprow, mypcol, m, n,  lmm,lnn,
                     isizeT, lrT1,lcT1,   lrB1,lcB1,   iia,jja,  iib,jjb );
               };
               CHKERR(cu_status);
             #endif
             };

          };
 

     }; /* end for jstart */

   }; 

  if (T != 0) {
    FreeHost( T );
    T = 0;
  };

  return;
}

#ifdef __cplusplus
extern "C"
#endif
void Cpdgecopy_d2h(int m, int n, 
        double *A, int ia, int ja, int *descA,
        double *B, int ib, int jb, int *descB  )
{
  /*
   * Transfer device to host
   */

   int ltrue = (1 == 1);
   int is_device_to_host = ltrue;
   Cpdgecopy_general(m,n,A,ia,ja,descA,B,ib,jb,descB, is_device_to_host );
} 

#ifdef __cplusplus
extern "C"
#endif
void Cpdgecopy_h2d(int m, int n, 
        double *A, int ia, int ja, int *descA,
        double *B, int ib, int jb, int *descB  )
{
  /*
   * Transfer host to device
   */
   int lfalse = (1 == 0);
   int is_device_to_host = lfalse;
   Cpdgecopy_general(m,n,A,ia,ja,descA,B,ib,jb,descB, is_device_to_host );
} 
