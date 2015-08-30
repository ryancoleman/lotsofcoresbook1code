#include "ooclu.h"

#ifdef __cplusplus
extern "C" 
#endif
int Cialast_simple( int ia, int mb, int myprow, int rsrc, int nprow )
{
  /*
   * return the last global index that belongs to "myprow"
   * but is less than or equal to "ia"
   */
  int iproc;
  int is_found;
  int iia;

  iia = -1;
  for(int i=0; i < nprow*mb; i++) {
    iia = ia - i;
    iproc = Cindxg2p( iia, mb,myprow,rsrc,nprow);
    is_found = (iproc == myprow);
    if (is_found) {
      break;
    };
  };

  return(iia);
}

#ifdef __cplusplus
extern "C" 
#endif
int Cialast_loop( int ia, int mb, int myprow, int rsrc, int nprow)
{

  int iia,ia_end; 

  iia = ia;
  while (myprow != Cindxg2p(iia,mb,myprow,rsrc,nprow) ) {
      ia_end = (iia-1) + (mb - MOD(iia-1,mb));
      iia = ia_end - mb;
  };


  if (iia >= 1) {
    assert( Cindxg2p(iia,mb,myprow,rsrc,nprow) == myprow );
  };
  if (iia-1 >= 1) {
    assert( Cindxg2p(iia-1,mb,myprow,rsrc,nprow) != myprow );
  };

  if (iia >= 1) {
     return( iia );
  }
  else {
    return(-1);
  };


}



#ifdef __cplusplus
extern "C" 
#endif
int Cialast_div( int ia, int mb, int myprow, int rsrc, int nprow )
{
  /*
   * return the last global index that belongs to "myprow"
   * but is less than "ia"
   */

  int iia, iproc;
  int is_found = 0;
  int myprow_shift, iproc_shift;

  iia = ia;
  iproc = Cindxg2p( iia,mb,myprow,rsrc,nprow);
  is_found = (iproc == myprow);
  if (is_found) { return(iia); };

  /*
   * shift rsrc
   */
  myprow_shift = myprow - rsrc;
  while (myprow_shift < 0) { myprow_shift += nprow; };
  iproc_shift = iproc - rsrc;
  while (iproc_shift < myprow) { iproc_shift += nprow; };


  iia =  ia -((iproc_shift-1)-myprow_shift)*mb ;
  iia = (iia-1) + (mb - MOD(iia-1,mb));

  /*
   * check
   */

  if (iia >= 1) {
     myprow = MOD(myprow,nprow);
     iproc = Cindxg2p( iia,mb,myprow,rsrc,nprow);
     assert( iproc == myprow );
     assert( Cindxg2p( iia+1,mb,myprow,rsrc,nprow) != myprow );
     assert( iia <= ia );

     return( iia );
     }
  else {
    return( -1 );
  };


}

#ifdef __cplusplus
extern "C" 
#endif
int Cialast( int ia, int mb, int myprow, int rsrc, int nprow )
{
  return( Cialast_simple(ia,mb,myprow,rsrc,nprow) );
}
