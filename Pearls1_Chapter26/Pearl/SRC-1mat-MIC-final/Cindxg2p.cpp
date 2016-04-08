#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
int Cindxg2p( int indxglob, int nb, int iproc, int isrcproc, int nprocs )
{
  int indxg2p;

  assert( nb >= 1);
  assert( nprocs >= 1);
  assert( (0 <= iproc) && (iproc < nprocs) );
  assert( (0 <= isrcproc) && (isrcproc < nprocs ) );

  indxg2p = ( MOD(isrcproc + (indxglob-1)/nb, nprocs ) );
  return( indxg2p );
}
