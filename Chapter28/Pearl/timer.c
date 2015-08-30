#include "timer.h"

#define CLK CLOCK_REALTIME
//#define CLK CLOCK_PROCESS_CPUTIME_ID
//#define CLK CLOCK_MONOTONIC

void timer_start(tv t) {
  clock_gettime(CLK,&t[0]);
}

void timer_stop(tv t) {
  clock_gettime(CLK,&t[1]);
  if(t[0].tv_nsec>t[1].tv_nsec) {
    t[1].tv_nsec+=1000000000;
    t[1].tv_sec--;
  }
  t[2].tv_sec+=(t[1].tv_sec-t[0].tv_sec);
  t[2].tv_nsec+=(t[1].tv_nsec-t[0].tv_nsec);
}

void timer_reset(tv t) {
  t[0].tv_sec=0;
  t[0].tv_nsec=0;
  t[1].tv_sec=0;
  t[1].tv_nsec=0;
  t[2].tv_sec=0;
  t[2].tv_nsec=0;
}

long int timer_nsec(tv t) {
  return 1000000000*t[2].tv_sec+t[2].tv_nsec;
}

double timer_sec(tv t) {
  return (double)timer_nsec(t)/1.0e+9;
}

double timer_msec(tv t) {
  return (double)timer_nsec(t)/1.0e+6;
}

double timer_usec(tv t) {
  return (double)timer_nsec(t)/1.0e+3;
}


#ifdef TIMER_TEST
main() {
  tv ts;
  timer_reset(ts);

  timer_start(ts);
  usleep(3000000);
  timer_stop(ts);

  printf(" sec: %lf\n",timer_sec(ts));
  printf("msec: %lf\n",timer_msec(ts));
  printf("usec: %lf\n",timer_usec(ts));
  printf("nsec: %ld\n",timer_nsec(ts));

  timer_start(ts);
  usleep(3000000);
  timer_stop(ts);

  printf(" sec: %lf\n",timer_sec(ts));
  printf("msec: %lf\n",timer_msec(ts));
  printf("usec: %lf\n",timer_usec(ts));
  printf("nsec: %ld\n",timer_nsec(ts));
}
#endif
