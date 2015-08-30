/*
 * common.h
 *
 *  Created on: Apr 30, 2013
 *      Author: lfeng
 */

#ifndef COMMON_H_
#define COMMON_H_

//#include "DebugLog.h"

#include <sys/time.h>
#include <iostream>

#undef __noinline
#undef __forceinline
#define __noinline             __attribute__((noinline))
#define __forceinline          inline __attribute__((always_inline))
#define __align(...)           __attribute__((aligned(__VA_ARGS__)))


#define TSLOG_DEBUG(message) std::cout << message << std::endl;
#define TSLOG_INFO(message) std::cout << message << std::endl;
#define TSLOG_ERROR(message) std::cout << message << std::endl;

__forceinline double dtime()
{

    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime, (struct timezone*) 0);
    tseconds = (double) (mytime.tv_sec +
                         mytime.tv_usec * 1.0e-6);
    return (tseconds);
}

#define CHECK_COI_RESULT(_COIFUNC) \
{ \
    COIRESULT result = _COIFUNC; \
    if (result != COI_SUCCESS) \
    { \
        printf("%s returned %s\n", #_COIFUNC, COIResultGetName(result));\
        return -1; \
    } \
}

#endif /* COMMON_H_ */
