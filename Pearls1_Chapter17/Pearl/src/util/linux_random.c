/*-----------------------------------------------------*\
 $Id: linux_random.c 20377 2011-05-24 13:29:18Z d3p852 $
\*-----------------------------------------------------*/
#include <stdlib.h>
#include "typesf2c.h"

void linux_sran_(Integer* input_seed)
{
  unsigned int seed;

  seed = (unsigned) *input_seed;
  (void) srandom(seed);
}
double linux_rand_(void)
{
  return (double) (((double) random())/(double) RAND_MAX);
}
