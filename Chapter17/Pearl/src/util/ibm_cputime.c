/*$Id: ibm_cputime.c 19707 2010-10-29 17:59:36Z d3y133 $*/

#if defined(USE_CLOCK)
#include <time.h>

double ibm_cputime_()
{
  return clock()/CLOCKS_PER_SEC;
}
#else
#include <sys/types.h>
#include <sys/times.h>
#define FAC2SEC 0.01

double ibm_cputime_()
{
  struct tms bufff;
  double e;
  time_t tt;
  double tarray[2];

  tt = times(&bufff);
  tarray[0] = FAC2SEC * bufff.tms_utime;
  tarray[1] = FAC2SEC * bufff.tms_stime;
  e = tarray[0] + tarray[1] ;
  return e;
}
#endif
