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

#include <mpi.h>
#include "ooclu.h"
#include <mkl_blacs.h>
//#include <mkl_scalapack.h>
//#include <mkl_pblas.h>

extern "C"
void pdmatgen(int *ICTXT, char *AFORM, char *DIAG,
              int *M, int *N, int *MB, int *NB, double *A,
              int *LDA, int *IAROW, int *IACOL, int *ISEED,
              int *IROFF, int *IRNUM, int *ICOFF, int *ICNUM,
              int *MYROW, int *MYCOL, int *NPROW, int *NPCOL);

extern "C"
void pdfillpad_(int *ICTXT, int *M, int *N, double *A, int *LDA,
               int *IPRE, int *IPOST, double *CHKVAL);

extern "C"
void pdchekpad_(int *ICTXT, char *MESS, int *M, int *N, double *A, int *LDA,
               int *IPRE, int *IPOST, double *CHKVAL);

extern "C"
void igsum2d_(int *ICTXT, char *SCOPE, char *TOP,
              int *M, int *N, int *A, int *LDA, int *RDEST, int *CDEST);

int pdlltinfo(int* n, int* nb, int* nprow, int* npcol, int* memsize);

int main(int argc, char** argv){
    int info; 
    int iam, nprocs, mydevice;
    int n, nb, nprow, npcol, memsize;
    int ICTXT, myprow, mypcol, lm, ln;
    int isroot;
    int i_one = 1, i_zero = 0, i_negone = -1;
    double d_one = 1.0, d_zero = 0.0, d_negone = -1.0;
    char *uplo = "Lower";
    int IASEED = 85;
    double *A;
    int descA[DLEN_];
    int prepad, midpad, postpad;
    double padval = -9923.0;
    double llttime;
    double mflops;
    size_t worksize;
    size_t sizebyte;

    blacs_pinfo_(&iam, &nprocs);
    isroot = (iam == 0);
    if(isroot) printf("read input\n"); ///////////////////////////////
    info = pdlltinfo(&n, &nb, &nprow, &npcol, &memsize);
    switch(info){
    case 1:
        if(isroot) printf("unable to read value of n\n");
        exit(0);
    case 2:
        if(isroot) printf("unable to read value of nb\n");
        exit(0);
    case 3:
        if(isroot) printf("unable to read value of p\n");
        exit(0);
    case 4:
        if(isroot) printf("unable to read value of q\n");
        exit(0);
    case 5:
        if(isroot) printf("unable to read value of memsize\n");
        exit(0);
    }
    if(nprow*npcol > nprocs){
        if(isroot) printf("grid size exceeded acceptable range\n");
        exit(0);
    }else if(nprow*npcol < nprocs){
        if(isroot) printf("too many processes\n");
        exit(0);
    }
    if(isroot){
        printf("the following values will be used\n");
        printf("n nb p q memsize\n");
        printf("%d %d %d %d %d\n",
            n, nb, nprow, npcol, memsize);
    }

    if(isroot) printf("grid init\n"); ////////////////////////////////
    blacs_get_(&i_negone, &i_zero, &ICTXT);
    blacs_gridinit_(&ICTXT, "R", &nprow, &npcol);
    blacs_gridinfo_(&ICTXT, &nprow, &npcol, &myprow, &mypcol);
    blacs_barrier_(&ICTXT, "All");

    if(isroot) printf("allocate memory\n"); //////////////////////////
    llttime = MPI_Wtime();

    lm = Cnumroc(n, nb, myprow, 0, nprow);
    ln = Cnumroc(n, nb, mypcol, 0, npcol);

//  prepad = MAX(lm, nb);
//  midpad = nb;
//  postpad = MAX(ln, nb);
    prepad = 0;
    midpad = 0;
    postpad = 0;

    Cdescinit(descA, n, n, nb, nb, i_zero, i_zero, ICTXT, lm+midpad, &info);

    // sizebyte = (prepad+(lm+midpad)*ln+postpad)*sizeof(double);
    worksize = lm+midpad;
    worksize = worksize*ln;
    worksize = worksize + prepad + postpad;
    sizebyte = worksize * sizeof(double);
    A = (double*) malloc(sizebyte);

    //check for error
    info = (A == NULL)? 1 : 0;
    igsum2d_(&ICTXT, "All", " ", &i_one, &i_one,
             &info, &i_one, &i_negone, &i_negone);
    if(info != 0){
        if(isroot) printf("unable to allocate\n");
        printf("process %d local work size %zd byte\n", iam, sizebyte);
        goto EXIT;
    }
    printf("process %d local work size %zd byte\n", iam, sizebyte);
/*  for(int i = 0; i < ln; i++){
        dscal(&lm, &d_zero, A+((size_t)i)*(size_t)lm, &i_one);
    }*/

/*  for(size_t df = 0; df < worksize; df++){
        A[df] = padval;
    }*/
    A = A + prepad;
    blacs_barrier_(&ICTXT, "All");
    if(isroot) printf("allocation time %lf\n", MPI_Wtime() - llttime);
    if(isroot) printf("generate A\n"); ///////////////////////////////
    llttime = MPI_Wtime();
        pdmatgen(&ICTXT, "Symm", "Diag", &descA[M_],
             &descA[N_], &descA[MB_], &descA[NB_],
             A, &descA[LLD_], &descA[RSRC_],
             &descA[CSRC_], &IASEED, &i_zero, &lm, &i_zero, &ln,
             &myprow, &mypcol, &nprow, &npcol);   
    blacs_barrier_(&ICTXT, "All");
    if(isroot) printf("generation time %lf\n", MPI_Wtime() - llttime);
/*
    if(isroot) printf("fill padding\n"); /////////////////////////////
    pdfillpad_(&ICTXT, &lm, &ln, A-prepad, &descA[LLD_],
               &prepad, &postpad, &padval);
    blacs_barrier_(&ICTXT, "All");
*/
/*  set 4*identity ///////////////////////////////////////////////////
    if(myprow==mypcol){
        for(int i=0; i<lm; i++){
            for(int j=0; j<i; j++){
                A[j*lm+i] = 0.0;
            }
            A[i*lm+i] = 4.0;
        }
    }else{
        for(int i=0; i<lm; i++){
            for(int j=0; j<ln; j++){
                A[j*lm+i] = 0.0;
            }
        }
    }
    printf("%d %d %lf %lf %lf\n", myprow, mypcol, A[0], A[lm+1], A[lm*ln-1]);
    blacs_barrier_(&ICTXT, "All");
*/

    #ifdef USE_MIC
        #ifdef __INTEL_OFFLOAD
            if(isroot)
                printf("offload compilation enabled\ninitialize each MIC\n");
            offload_init(&iam, &mydevice);
        #else
            if(isroot)
                printf("offload compilation not enabled\n");
            exit(0);
        #endif
    #else
        cublasInit();
    #endif
    profinit();
    if(isroot) printf("LLT\n"); //////////////////////////////////////
    llttime = MPI_Wtime();
    pdpotrf_ooc2(uplo,&n,A,&i_one,&i_one,descA,&memsize,&info);
    //scalapack_pdpotrf(uplo,&n,A,&i_one,&i_one,descA,&info);
    blacs_barrier_(&ICTXT, "All");
    llttime = MPI_Wtime() - llttime;
    if(isroot) profstat();
    blacs_barrier_(&ICTXT, "All");
    if(myprow+mypcol == 1) profstat();
    blacs_barrier_(&ICTXT, "All");
    if(isroot){
        mflops = (double) n;
        mflops = mflops/3.0 + 0.5;
        mflops = mflops/llttime/1024.0/1024.0;
        mflops = mflops*(double)(n)*(double)(n);
        printf("result: n nb p q memsize(K wordsize) llt-time MFLOPS\n");
        printf("test result: MIC %d %d %d %d %d %lf %lf\n",
            n, nb, nprow, npcol, memsize/1024, llttime, mflops);
    }
    blacs_barrier_(&ICTXT, "All");
/*
    if(isroot) printf("check padding\n"); ////////////////////////////
    pdchekpad_(&ICTXT, "LLT", &lm, &ln, A-prepad, &descA[LLD_],
              &prepad, &postpad, &padval);
    blacs_barrier_(&ICTXT, "All");
*/
/*
    printf("%d %d %lf %lf %lf\n", myprow, mypcol, A[0], A[lm+1], A[lm*ln-1]);
    blacs_barrier_(&ICTXT, "All");
*/

    free(A-prepad);
    if(isroot) printf("driver complete\n"); //////////////////////////
    blacs_gridexit_(&ICTXT);
    #ifdef USE_MIC
        offload_destroy();
    #else
        cublasShutdown();
    #endif
    EXIT:
    blacs_exit_(&i_zero);
    return 0;
}
