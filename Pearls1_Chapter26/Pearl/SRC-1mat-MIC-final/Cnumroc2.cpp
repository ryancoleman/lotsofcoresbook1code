#include "ooclu.h"

#ifdef __cplusplus
extern "C"
#endif
int Cnumroc2(int ia,int n, int nb, int iproc, int isrcproc,int nprocs )
{

int m1, m2;


m1 = Cnumroc( ia, nb,iproc,isrcproc,nprocs );
m2 = Cnumroc( ia+n-1, nb,iproc,isrcproc,nprocs );

return( m2-m1+1 );
}
