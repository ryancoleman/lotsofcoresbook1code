#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
void Cinfog2l( int grindx, int gcindx, int *desc, 
         int nprow,int npcol, int myrow, int mycol,
         int *lrindx, int *lcindx, int *rsrc, int *csrc )
{
  scalapack_infog2l( &grindx, &gcindx, desc,
               &nprow, &npcol, &myrow, &mycol,
               lrindx, lcindx, rsrc, csrc );
}
