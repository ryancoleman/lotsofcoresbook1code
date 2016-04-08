#include "ooclu.h"
#ifdef __cplusplus
extern "C"
#endif


#include "mpi.h"

void pdcopymatrix(double *  	A, int   	ia, int   	ja, int *  	descA)
{
    int myid, numproc, i, icontxt;;
    int Locp, Locq, lld, info;
    int mb,nb, ic_proc, jc_proc;
    int ctxt, nprow,npcol,myprow,mypcol;
    char scope,top;
    int ldAtmp;
    int isizedAtmp;
    double *dAtmp, *hAtmp, *hcopydAtmp;
    int cuerr,hosterr;
    const int elementSize = sizeof(double);
    const int nelements = 100;

    MPI_Comm comm;
       
    comm=MPI_COMM_WORLD;
    icontxt = descA[CTXT_];
    Cblacs_gridinfo( icontxt, &nprow,&npcol,&myprow,&mypcol);
    assert( nprow >= 1);
    assert( npcol >= 1);
    assert( (0 <= myprow) && (myprow < nprow));
    assert( (0 <= mypcol) && (mypcol < npcol));

/**************  Set up device matrix *****************/
    jc_proc = 0;
    isizedAtmp = nelements ;
    isizedAtmp *= elementSize;

    cuerr = cudaMalloc((void**) &dAtmp, isizedAtmp);
    hAtmp = (double*)malloc(isizedAtmp);
    assert( dAtmp !=0 );

/**  assign value to the matrix on device */
    if ( mypcol == jc_proc ) {
      for (i=0;i<nelements;i++) {
        hAtmp[i]=1.02;
      }
    } else {
      for (i=0;i<nelements;i++) {
        hAtmp[i]=2.0;
      }
    }
    for (i=0;i<nelements;i++) {
        hcopydAtmp[i]=0.001;
    }
    cudaMemcpy (dAtmp, hAtmp, isizedAtmp, cudaMemcpyHostToDevice);
    free (hAtmp);

//    printf("after cudaMemset,isizedAtmp=%d\n",isizedAtmp);

/*****************Set up device matrix *****************/          
            scope = 'R';
            top =' ';
//           Locp = isizedAtmp;
            Locp = nelements;
            Locq = 1;
            lld = Locp;

            printf("before dgebs/r2d: icontxt=%d,Locp=%d,Locq=%d,lld=%d,\n",
            icontxt,Locp,Locq,lld);

//            MPI_Bcast(dAtmp, nelements, MPI_DOUBLE, 0, comm); 


            if (mypcol == jc_proc) {
                scalapack_dgebs2d(&icontxt, &scope, &top,
                                  &Locp, &Locq, dAtmp, &lld );
            } else {
                scalapack_dgebr2d(&icontxt,&scope,&top,
                                  &Locp,&Locq, dAtmp, &lld,
                                  &myprow, &jc_proc );
            };


     cudaMemcpy (hcopydAtmp, dAtmp, isizedAtmp, cudaMemcpyDeviceToHost);
     if (mypcol == 1 && myprow == 1) {
     printf( "mypcol %d, myprow %d \n", mypcol, myprow);
            for (int i=0; i<nelements; i=i+4){
              printf( "%f, %f, %f, %f \n", 
//                      dAtmp[i],dAtmp[i+1],dAtmp[i+2],dAtmp[i+3]);
                      hcopydAtmp[i],hcopydAtmp[i+1],hcopydAtmp[i+2],hcopydAtmp[i+3]);
            }
     }
    cudaFree(dAtmp);
}

#ifdef __cplusplus
extern "C"
#endif
void pdcopymatrix_(double *A, int ia, int ja, int *descA)
{
  pdcopymatrix(A, ia, ja, descA);
}
