#ifndef SCALAPACK_I_H
#define  SCALAPACK_I_H 1

#if defined(ADD_)
#define scalapack_picopy  picopy_
#define scalapack_pilaprnt pilaprnt_

#define scalapack_numroc  numroc_
#define scalapack_indxg2p  indxg2p_
#define scalapack_infog2l  infog2l_
#define scalapack_descinit  descinit_

#define scalapack_igebr2d igebr2d_
#define scalapack_igebs2d igebs2d_
#define scalapack_igsum2d igsum2d_
#define scalapack_pielget pielget_
#else

#define scalapack_picopy  picopy
#define scalapack_pilaprnt pilaprnt

#define scalapack_numroc  numroc
#define scalapack_indxg2p  indxg2p
#define scalapack_infog2l  infog2l
#define scalapack_descinit  descinit

#define scalapack_igebr2d igebr2d
#define scalapack_igebs2d igebs2d
#define scalapack_igsum2d igsum2d
#define scalapack_pielget pielget

#endif

extern "C" {

  int scalapack_numroc( int *n, int *nb, 
                  int *iproc, int *isrcproc, int *nprocs );

  int scalapack_indxg2p( int *indxglob, int *nb, 
                   int *iproc, int *isrcproc, int *nprocs);

  void scalapack_infog2l( int *grindx, int *gcindx, int *desc, 
                 int *nprow, int *npcol, int *myrow, int *mycol,
                 int *lrindx, int *lcindx, int *rsrc, int *csrc );

  void scalapack_descinit( int *desc, int *m, int *n, int *mb, int *nb,
                  int *irsrc, int *icsrc, int *ctxt, int *lld, int *info );



void Cblacs_gridinfo( int icontxt, int *nprow, int *npcol, 
                                   int *myprow, int *mypcol);



void scalapack_picopy( int *n, 
           int *X, int *ix, int *jx, int *descX, int *incX,
           int *Y, int *iy, int *jy, int *descY, int *incY );


void scalapack_pilaprnt( int *m, int *n, 
       int *A, int *ia, int *ja, int *descA, 
       int *irprnt, int *icprnt, char *cmatnm, int *nout, int *work );


void scalapack_igebs2d( int *icontxt, char *scope, char *top,
           int *m, int *n, int *A, int *lda );

void scalapack_igebr2d( int *icontxt, char *scope, char *top,
           int *m, int *n, int *A, int *lda,
           int *rsrc, int *csrc );



void scalapack_igsum2d( int *icontxt, char *scope, char *top,
       int *m, int *n, int *A, int *lda,    int *rdest, int *cdest );



void scalapack_pielget( char *scope, char *top, int *alpha, 
    int *A, int *ia, int *ja, int *descA );

}
#endif
