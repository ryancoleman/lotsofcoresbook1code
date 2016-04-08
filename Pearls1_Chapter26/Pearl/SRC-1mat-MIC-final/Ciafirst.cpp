#include "ooclu.h"

#ifdef __cplusplus
extern "C" 
#endif
int Ciafirst_simple(int ia, int mb, int myprow, int rsrc, int nprow) 
{
  int iia,iproc;
  int is_found;

  for(iia=ia; iia <= ia + mb*nprow; iia++) {
    iproc = Cindxg2p(iia,mb,myprow,rsrc,nprow);
    is_found = (iproc == myprow);
    if (is_found) {
      break;
    };
  };

  assert( Cindxg2p(iia,mb,myprow,rsrc,nprow) == myprow );
  return(iia);
}


#ifdef __cplusplus
extern "C" 
#endif
int Ciafirst_loop(int ia, int mb, int myprow, int rsrc, int nprow) 
{

  const int idebug = 1;
  int iproc;
  int ia_first, ia_end;
  int is_found;


  ia_first = ia;
  iproc = Cindxg2p(ia_first,mb,myprow,rsrc,nprow);
  is_found = (iproc == myprow);
  if (is_found) {
    return( ia_first );
  };

  for(int i=1;  (!is_found) && (i <= nprow+1); i++)  {
    iproc = Cindxg2p(ia_first,mb,myprow,rsrc,nprow);
    is_found = (iproc == myprow);
    if (is_found) { 
      break;
    };

    ia_end =  (ia_first-1) + (mb - MOD( ia_first-1, mb ));
    assert( MOD(ia_end,mb) == 0);

    ia_first = ia_end + 1;
  };

  if (idebug >= 1) {
    assert( Cindxg2p(ia_first,mb,myprow,rsrc,nprow) == myprow);
    if (ia_first > ia) {
      assert( Cindxg2p(ia_first-1,mb,myprow,rsrc,nprow) != myprow);
    };
  };

  return( ia_first );
}




#ifdef __cplusplus
extern "C" 
#endif
int Ciafirst_div( int ia, int mb, int myprow, int rsrc, int nprow )
{
  /*
   * return the 1st global index that belongs to "myprow"
   * that is higher than "ia"
   */
  const int idebug = 1;

  int iia,iproc,iproc0;
  int is_found = 0;
  int myprow_shift, iproc_shift;

  iia = ia;
  iproc0 = Cindxg2p( iia,mb,myprow,rsrc,nprow);
  iproc = iproc0;
  is_found = (iproc0 == myprow);
  if (is_found) { return( iia ); };


  iproc_shift = iproc - rsrc;
  while (iproc_shift < 0) { iproc_shift += nprow; };

  myprow_shift = myprow - rsrc;
  while (myprow_shift < iproc_shift) { myprow_shift += nprow; };


  iia = ia + (mb - MOD( (ia-1),mb)) + ((myprow_shift-1)-iproc_shift)*mb;

  /*
   * Check result
   */
  if (idebug >= 1) {

  iproc = Cindxg2p( iia,mb,myprow,rsrc,nprow );
  assert( iproc == myprow );

  assert( iia > ia );
  if (MOD(iia-1,mb)!=0) {
    printf("ia %d iia %d mb %d rsrc %d myprow %d nprow %d\n",
            ia,   iia,   mb,   rsrc,   myprow,   nprow );
  };
  assert( MOD(iia-1,mb) == 0);
  assert( Cindxg2p( iia-1,mb,myprow,rsrc,nprow) != myprow );
  };


  return( iia );
}

#ifdef __cplusplus
extern "C" 
#endif
int Ciafirst(int ia, int mb, int myprow, int rsrc, int nprow) 
{
  return(Ciafirst_div(ia,mb,myprow,rsrc,nprow));
}
