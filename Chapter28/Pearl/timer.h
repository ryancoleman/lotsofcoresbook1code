#ifndef TIMER_H
#define TIMER_H

#include <time.h>

typedef struct timespec tv[3];

void timer_start(tv t);
void timer_stop(tv t);
void timer_reset(tv t);
long int timer_nsec(tv t);
double timer_sec(tv t);
double timer_msec(tv t);
double timer_usec(tv t);

#endif
