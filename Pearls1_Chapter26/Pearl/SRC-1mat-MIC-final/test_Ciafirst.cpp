#include "ooclu.h"

#ifdef __cplusplus
extern "C" {
#endif
  int Ciafirst_simple(int ia, int mb, int myprow, int rsrc, int nprow );
  int Ciafirst_loop(int ia, int mb, int myprow, int rsrc, int nprow );
  int Ciafirst_div(int ia, int mb, int myprow, int rsrc, int nprow );
#ifdef __cplusplus
}
#endif


main()
{

  int nprow,rsrc,mb,ia,myprow;
  int iproc_simple, iproc_loop, iproc_div;
  int nerror_loop;
  int nerror_div;

  int mb_max, nprow_max;


  nerror_loop = 0;
  nerror_div = 0;

  nprow_max = 15;
  mb_max = 128;


  for(nprow=0; nprow <= nprow_max; nprow++) {
  for(myprow=0; myprow <= (nprow-1); myprow++) {
  for(rsrc=0; rsrc <= (nprow-1); rsrc++) {
  for(mb=mb_max; mb <= mb_max; mb++) {
  for(ia=1; ia <= 2*nprow*mb; ia++) {
        iproc_simple = Ciafirst_simple(ia,mb,myprow,rsrc,nprow);
        iproc_loop = Ciafirst_loop(ia,mb,myprow,rsrc,nprow);
        iproc_div = Ciafirst_div(ia,mb,myprow,rsrc,nprow);

        if (iproc_simple != iproc_loop) {
            nerror_loop++;
        };
        if (iproc_simple != iproc_div) {
            nerror_div++;
        };
  };
  };
  };
  };
  };

  printf("nerror_loop %d nerror_div %d \n",
          nerror_loop,   nerror_div );

}


