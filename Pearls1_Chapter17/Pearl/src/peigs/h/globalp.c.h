/*
 $Id: globalp.c.h 19708 2010-10-29 18:04:21Z d3y133 $

 file globalp.c.h 
 */

#ifdef  STD_INT
typedef int    Integer;
typedef unsigned int unInteger;
#else
typedef long   Integer;
typedef unsigned long unInteger;
#endif

#ifdef  STD_DBL
typedef   double         DoublePrecision;
#else
typedef   long double    DoublePrecision;
#endif

#include "blas_lapack.h"
#include "peigs_types.h"

