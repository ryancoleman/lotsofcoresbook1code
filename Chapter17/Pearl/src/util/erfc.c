#include <math.h>

extern double erfc (double);

double derfc_ (double * x)
{
  return (erfc (*x));
}

double erfc_ (double * x)
{
  return (erfc (*x));
}
/* $Id: erfc.c 21176 2011-10-10 06:35:49Z d3y133 $ */
