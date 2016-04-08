#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
int Cnumroc( int n, int nb, int iproc, int isrcproc, int nprocs )
{
int mydist, nblocks, numroc, extrablks;

mydist = MOD( nprocs+iproc-isrcproc,nprocs);
nblocks = n / nb;
numroc = (nblocks/nprocs) * nb;
extrablks = MOD( nblocks, nprocs );

if (mydist < extrablks) {
    numroc = numroc + nb;
    }
else if (mydist == extrablks) {
   numroc = numroc + MOD( n,nb);
   };

return( numroc );
}
