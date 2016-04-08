#if defined(BL_FORT_USE_UNDERSCORE)
#define WALL_SECOND wall_second_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define WALL_SECOND wall_second__
#elif defined(BL_FORT_USE_LOWERCASE)
#define WALL_SECOND wall_second
#endif

#if defined(_BL_ANSI_TIME)

#include <stdlib.h>
#include <time.h>

void
WALL_SECOND(double* rslt)
{
  static time_t start;
  static int inited = 0;
  if ( inited == 0 )
    {
      inited = 1;
      start = time(0);
    }
  *rslt = (double)(time(0)-start);
}

#elif defined(WIN32)

#include <windows.h>
#include <time.h>

void
WALL_SECOND(double* rslt)
{
  static int inited = 0;
  static LONGLONG llStart;
  static double rate;
  LARGE_INTEGER li;

  if ( inited == 0 )
    {
      LARGE_INTEGER lii;
      inited = 1;
      QueryPerformanceFrequency(&lii);
      rate = 1.0/lii.QuadPart;
      QueryPerformanceCounter(&lii);
      llStart = lii.QuadPart;
    }
  QueryPerformanceCounter(&li);
  *rslt = (double)((li.QuadPart-llStart)*rate);
}

#else

#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if defined(_BL_USE_MPI_WTIME)
#include <mpi.h>

void
WALL_SECOND(double* rslt)
{
  int fr;
  assert(  MPI_Initialized(&fr) == 0 );
  *rslt = MPI_Wtime();
}
#else
void
WALL_SECOND(double* rslt)
{
  struct timeval tp;
  /* cannot fail, so why check */
  gettimeofday(&tp, 0);
  *rslt = (double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6;
}
#endif
#endif

#include <stdlib.h>

#if defined(BL_FORT_USE_UNDERSCORE)
#define SYS_ABORT sys_abort_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define SYS_ABORT sys_abort__
#elif defined(BL_FORT_USE_LOWERCASE)
#define SYS_ABORT sys_abort
#endif

void
SYS_ABORT()
{
  abort();
}
