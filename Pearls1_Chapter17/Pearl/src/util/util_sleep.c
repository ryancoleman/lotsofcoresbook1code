/*
 $Id: util_sleep.c 26140 2014-09-03 00:19:44Z edo $
 */

#if defined(CRAY) &&!defined(__crayx1)

void UTIL_SLEEP(long *t)
{}

#elif defined(WIN32) && !defined(__MINGW32__)

#include "winutil.h"
#include "typesf2c.h"
void FATR UTIL_SLEEP(long *t) 
{
  unsigned s = *t;
  sleep(s);
}

#else

#include <unistd.h>
#if defined(__MINGW32__)
#  include <windows.h>
#  define sleep(x) Sleep(1000*(x))
#endif
void util_sleep_(long *t) 
{
  unsigned s = *t;
  sleep(s);
}

#endif
