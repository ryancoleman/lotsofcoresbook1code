#include "ooclu.h"

void Cpzlaprnt( int m, int n, double *A, int ia, int ja, int *descA, char *cmatnm )
{
  double zalpha[] =  {0.0, 0.0};

  char scope = 'A';
  char top = ' ';
  int i,j,iia,jja;
  int nprow,npcol,myprow,mypcol;
  int is_root;
  int is_ok;


  Cblacs_gridinfo( descA[CTXT_], &nprow,&npcol, &myprow,&mypcol);
  is_root = (myprow == 0) && (mypcol == 0);
  if (is_root) {
    printf("Cpzlaprnt: m %d n %d ia %d ja %d cmatnm %s\n",
                       m,   n,   ia,   ja,   cmatnm );
  };

  for(i=1; i <= m; i++) {
    for(j=1; j <= n; j++) {
      iia = (ia-1) + i;
      jja = (ja-1) + j;

      is_ok = ( (1 <= iia) && (iia <= descA[M_]) ) &&
              ( (1 <= jja) && (jja <= descA[N_]) );
      if (is_root && (!is_ok)) {
        printf("i %d iia %d descA[M_] %d j %d jja %d descA[N_] %d \n",
                i,   iia,   descA[M_],   j,   jja,   descA[N_] );
      };
      assert( is_ok );

      scalapack_pzelget( &scope, &top, zalpha, A, &iia, &jja, descA );
      if (is_root) {
        printf("%s(%d,%d) = %lf + I * %lf\n",
               cmatnm, i,j,  zalpha[REAL_PART], zalpha[IMAG_PART] );
      };
    };
  };

}

