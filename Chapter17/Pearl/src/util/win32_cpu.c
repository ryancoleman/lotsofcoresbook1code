/*
 $Id: win32_cpu.c 19707 2010-10-29 17:59:36Z d3y133 $
*/

#include <time.h>
#include "typesf2c.h"
#ifdef WIN32
double FATR WIN32_CPUTIME(void)
#else
double win32_cputime_(void)
#endif
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
