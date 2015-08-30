*
* $Id: datasizes.h 19708 2010-10-29 18:04:21Z d3y133 $
*
*     Number of bytes in an INTEGER variable.
      INTEGER    NBYTEI
#ifdef STD_INT
      PARAMETER (NBYTEI = 4)
#else
      PARAMETER (NBYTEI = 8)
#endif
