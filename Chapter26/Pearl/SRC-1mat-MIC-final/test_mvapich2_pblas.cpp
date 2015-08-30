#include "ooclu.h"
//#ifdef __cplusplus
extern "C"
//#endif
int main()
{
    int myid, numproc, i;
    int Locp, Locq, lld, info;
    int mb,nb, ic_proc, jc_proc;
    int ctxt, nprow,npcol,myprow,mypcol;
    char scope,top;
    int ldAtmp;
    int isizedAtmp;
    double *dAtmp;
    int cuerr;
    const int elementSize = sizeof(double);
    const int nelements = 100;
/*
*        Define process grid
*/
    Cblacs_pinfo(&myid, &numproc);
    Cblacs_get(0, 0, &ctxt);
/*
  setup P x 2  grid
*/
    npcol =2; 
    nprow = numproc/npcol;
    Cblacs_gridinit(&ctxt, "Row-major", nprow, npcol);
    Cblacs_gridinfo( ctxt, &nprow,&npcol,&myprow,&mypcol);
    assert( nprow >= 1);
    assert( npcol >= 1);
    assert( (0 <= myprow) && (myprow < nprow));
    assert( (0 <= mypcol) && (mypcol < npcol));

/**************  Set up device matrix *****************/
    jc_proc = 0;
    isizedAtmp = nelements ;
    isizedAtmp *= elementSize;

    cuerr = cudaMalloc((void**) &dAtmp, isizedAtmp);
    assert( dAtmp !=0 );
    if ( mypcol == jc_proc ) {
            cudaMemset(dAtmp, 0, isizedAtmp);
    } else {
            cudaMemset(dAtmp, 1, isizedAtmp);
    }
/*****************Set up device matrix *****************/          
            scope = 'R';
            top =' ';
            Locp = isizedAtmp;
            Locq = 1;
            lld = Locp;
            if (mypcol == jc_proc) {
                scalapack_dgebs2d(&ctxt, &scope, &top,
                                  &Locp, &Locq, dAtmp, &lld );
            } else {
                scalapack_dgebr2d(&ctxt,&scope,&top,
                                  &Locp,&Locq, dAtmp, &lld,
                                  &myprow, &jc_proc );
            };
            fprintf(stdout, "myypcol %d, myprow %d \n", mypcol, myprow);
            for (int i=0; i<=nelements; i=i+4){
              fprintf(stdout, "%f, %f, %f, %f \n", 
                      dAtmp[i],dAtmp[i+1],dAtmp[i+2],dAtmp[i+3]);
            }
            fflush(stdout);
    cudaFree(dAtmp);
    return 0;
}

