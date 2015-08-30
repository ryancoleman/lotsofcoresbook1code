
#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
int Cindxl2g(int indxloc, int nb, int iproc, int isrcproc, int nprocs ) 
{
  int indxl2g;

  assert( nb >= 1);
  assert( nprocs >= 1);
  assert( (0 <= iproc) && (iproc < nprocs) );
  assert( (0 <= isrcproc) && (isrcproc < nprocs ) );

  indxl2g = nprocs*nb*((indxloc-1)/nb) + MOD(indxloc-1,nb) +
                          MOD(nprocs+iproc-isrcproc, nprocs)*nb + 1;
  return(indxl2g);


}
